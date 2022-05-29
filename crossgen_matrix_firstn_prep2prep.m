%% Cross-validation matrix of preparatory attention. Done separately within each ROI
% written by Tanya Wen
% 2019/01/23

dbstop if error

% addpath(genpath('/imaging/tw05'))
addpath('/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp')
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
%addpath(genpath('/imaging/local/software/spm_toolbox/eeglab13_4_3b'))
% if ~strcmp('EEG',spm_get_defaults('modality'))
spm('defaults', 'eeg');
% end
tempdir='/imaging/dm01/temp/'; % djm, temp location for me to write to

% required software: libSVM & RSAtoolbox
workingdir = '/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp';
addpath(genpath(fullfile('/imaging/tw05/Preparatory_Attention_Study/','software','rsatoolbox')));
addpath(fullfile('/imaging/tw05/Preparatory_Attention_Study/','software','libsvm-mat-2.87-1'));
% control variables
libSVMsettings='-s 1 -t 0'; % nu-SVM, linear
nRandomisations=1000;
nfolds = 5;
% rmpath('/hpc-software/matlab/r2014a/toolbox/bioinfo/biolearning/'); % to make sure libSVM code is used (not strictly necessary: matlab svmtrain yields exactly same model)

% Define SUBJECT INFORMATION
subs = [1,2,3,4,5,6,7,8,9,10,11,13,15,16,17,18,19,20];  % subject numbers
subjects_dirs = {'meg16_0317/161107','meg16_0319/161110','meg16_0321/161111','meg16_0322/161114','meg16_0325/161115','meg16_0327/161117','meg16_0330/161121','meg16_0332/161122','meg16_0333/161124','meg16_0337/161128','meg16_0339/161129','meg16_0340/161129','meg16_0341/161201','meg16_0343/161202','meg16_0345/161206','meg16_0346/161206','meg16_0348/161208','meg16_0349/161208','meg16_0350/161212','meg16_0352/161213'};
subjnum = [1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,6,2]; % counterbalancing numbers

min_trials_mvpa = [88,98,87,98,67,67,91,89,56,82,59,0,66,0,90,84,45,68,97,93];

%% parallelize
nw=100;
scheduler=cbu_scheduler();
scheduler.SubmitArguments='-q compute -l mem=250gb -l walltime=460800';

if isempty(gcp('nocreate')) || ~exist('pool','var') || pool.NumWorkers ~= nw,
    if ~isempty(gcp('nocreate'))
        delete(gcp('nocreate'))
    end
    scheduler.NumWorkers=nw;
    pool=parpool(scheduler,nw);
end

%%%%%%% djm, moved first bit of ROI loop outside subject loop
roi_folder = '/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp/ROIs';
roi_list = {'MD_ESV.nii'}; %{'Auditory_Te3.nii','Visual_hOc1.nii','LPFC.nii'};
roinames=roi_list;
roinames = strrep(roinames, '_', '-');
roinames = strrep(roinames, '.nii', '');

coords=cell(1,numel(roi_list));
for rois = 1:numel(roi_list)
    
    Vroi = spm_vol(fullfile(roi_folder,roi_list{rois}));
    [Y, XYZ]=spm_read_vols(Vroi);
    coords{rois} = XYZ(:,Y>0)';
end
%%%%%%%%

for sub = [2,3,4,5,6,7,8,9,10,11,13,15,16,17,18,19,20]
    
    tic
    fprintf('\nLoading data from subject %d...',sub)
    
    % go to subject directory
    cd(workingdir)
    swd = sprintf('sub%02d/%s',sub,subjects_dirs{sub}); % subject working directory
    cd(swd)
    
    
    %% attentional template
    % load data
    D = spm_eeg_load('caefMspm12_attention_task_block1_raw.mat');
    val = 1; %inversion method = MMN
    %vertex  = D.inv{val}.mesh.tess_mni.vert(D.inv{val}.inverse.Is, :); % find existing vertices
    
    %% define labels
    
    % define subject specific pairings
    allperms = perms([1 2 3]);
    for i = 1:length(allperms)
        if rem(subjnum(sub),6)+1 == i
            choosetargets = allperms(i,:);
        end
    end
    
    fprintf('Took %.1f minutes.',toc/60);
    
    for rois = 1:numel(roi_list)
        
        tic
        fprintf('\nExtracting data from ROI %d of subject %d...',rois,sub)
        
        % extract voxel x time x trials matrix of ROI
        D.val=val;
        D.inv{val}.source.XYZ = coords{rois};
        D.inv{val}.source.rad = 0;
        D.inv{val}.source.label= roinames(rois);
        D.inv{val}.source.fname = strcat(roinames{rois},'_attn1');
        D.inv{val}.source.type  = 'trials';
        try
            Dattn=spm_eeg_inv_extract_tw(D); % return header to avoid needing to reload
        catch
            D.inv{val}.source.fname = fullfile(tempdir,D.inv{val}.source.fname); % djm, in case don't have write permission to original directory
            Dattn=spm_eeg_inv_extract_tw(D);
        end
        
        data3D = single(Dattn.fttimelock.trial);
        cond_values = Dattn.fttimelock.trialinfo;
        
        % define cue1 & cue2 index vectors
        cues=cell(1,2);
        for i=1:2
            cues{i} = find(cond_values==i);
        end
        
        fprintf('Took %.1f minutes.',toc/60);
        tic
        fprintf('\nAnalysing ROI %d of subject %d',rois,sub)
        
        %% smoothing option
        smoothing = {'Raw','Smooth'};
        for smooth = 1%:numel(smoothing)
            
            %% sliding time window
            winsizes = 8;%[8,25,50,125]; %32ms,100ms,200ms,500ms
            for winsize = winsizes
                
                % set analysis time windows
                twin_train = [1 size(data3D,3)]; % time of interest
                stpsize        = 1; %step size
                trainwins      = [twin_train(1):stpsize:(twin_train(2)-winsize)]';
                trainwins(:,2) = trainwins(:,1) + winsize;
                ntrains        = size(trainwins,1);
                twin_test = [1 size(data3D,3)]; % time of interest
                testwins      = [twin_test(1):stpsize:(twin_test(2)-winsize)]'; %#ok<*NBRAK>
                testwins(:,2) = testwins(:,1) + winsize;
                ntests         = size(testwins,1);
                %%%%% djm: this bit allows the variable to be 'sliced' in parfor loop;
                windatatrain=nan(size(data3D,1),size(data3D,2),winsize+1,ntrains);
                % the window seems to be 8+1 samples long
                for islide = 1:ntrains
                    windatatrain(:,:,:,islide)=data3D(:,:,trainwins(islide,1):trainwins(islide,2));
                end
                windatatest=windatatrain; % making a separate variable stops it needing to be 'broadcast'
                %%%%
                
                accuracy_matrix = nan(ntrains,ntests);
                predicted_vals = cell(ntrains,ntests);
                
                
                
                tempind = [cues{:}];
                train_cond_values = cond_values(tempind);
                c1 = tempind(train_cond_values==1);
                c2 = tempind(train_cond_values==2);
                subsetind = [tempind(c1(1:min_trials_mvpa(sub)/2)),tempind(c2(1:min_trials_mvpa(sub)/2))];
                labels_train=cond_values(subsetind)';
                
                testind = setdiff([cues{:}],subsetind);
                labels_test=cond_values(testind)';
                
                parfor islide = 1:ntrains
                    patternsTrain=double(windatatrain(subsetind,:,:,islide));
                    if smooth == 1 % no smoothing
                        patternsTrain=reshape(patternsTrain,length(subsetind),[]);
                    elseif smooth == 2 % smoothing
                        patternsTrain = mean(patternsTrain,3);
                    end
                    
                    % classification using linear SVM (train the classifier)
                    model=svmtrain(labels_train,patternsTrain,libSVMsettings);
                    
                    for jslide = 1:ntests
                        itest = data3D(:,:,testwins(jslide,1):testwins(jslide,2));
                        patternsTest = double(itest(testind,:,:));
%                         patternsTest = double(windatatest(testind,:,:,jslide));
                        if smooth == 1 % no smoothing
                            patternsTest=reshape(patternsTest,length(testind),[]);
                        elseif smooth == 2 % smoothing
                            patternsTest = mean(patternsTest,3);
                        end
                        
                        % classification using linear SVM (test the classifier)
                        [predicted,accuracy,~]=svmpredict(labels_test,patternsTest,model);
                        accuracy_matrix(islide,jslide)=accuracy(1);
                        predicted_vals{islide,jslide} = predicted;
                        
                        
                    end % next slide
                end % next slide
                
                
                output_dir = 'crossgen_ROIs_subsample_prep-prep';
                try
                    cd(output_dir)
                catch
                    eval(sprintf('!mkdir %s',output_dir)); cd(output_dir);
                end
                
                figure(sub)
                imagesc(accuracy_matrix);
                colorbar;
                caxis([1,100]);
                
                % djm:
                set(gca,'xticklabel',(cellfun(@str2double,get(gca,'xticklabel'))-1)*stpsize*4-100);
                set(gca,'yticklabel',(cellfun(@str2double,get(gca,'yticklabel'))-1)*stpsize*4-100);
                xlabel('start of test windows (ms)')
                ylabel('start of train windows (ms)')
                
                fprintf('\nTook %.1f hours.',toc/60/60);
                
                saveas(gcf,sprintf('CrossGenMatrix_firstn_%03dms%s_%s_sub%s.png',winsize*4,smoothing{smooth},roinames{rois},num2str(sub)));
                save(sprintf('CrossGenMatrix_firstn_%03dms%s_%s_sub%s.mat',winsize*4,smoothing{smooth},roinames{rois},num2str(sub)));
                close;
                cd(workingdir)
                cd(swd)
                
            end % winsize
        end % smooth
    end % rois
end

delete(gcp)