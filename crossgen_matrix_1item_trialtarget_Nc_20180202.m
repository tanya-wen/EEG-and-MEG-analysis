cl%% Cross-validation of of which target the participant is attending to in 3-item displays. Done separately within each ROI
% written by Tanya Wen
% 2016/10/25
dbstop if error

addpath(genpath('/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp'))
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'))
addpath(genpath('/imaging/local/software/spm_toolbox/eeglab13_4_3b'))
spm('defaults', 'eeg');

% required software: libSVM & RSAtoolbox
workingdir = '/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp';
addpath(genpath(fullfile(workingdir,'software','rsatoolbox')));
addpath(fullfile('/imaging/tw05/Preparatory_Attention_Study','software','libsvm-mat-2.87-1'));
% control variables
libSVMsettings='-s 1 -t 0'; % nu-SVM, linear
nRandomisations=1000;
rmpath('/hpc-software/matlab/r2015a/toolbox/bioinfo/biolearning/'); % to make sure libSVM code is used (not strictly necessary: matlab svmtrain yields exactly same model)
nfolds = 5;

% Define SUBJECT INFORMATION
subs = [1,2,3,4,5,6,7,8,9,10,11,13,15,16,17,18,19,20];  % subject numbers
subjects_dirs = {'meg16_0317/161107','meg16_0319/161110','meg16_0321/161111','meg16_0322/161114','meg16_0325/161115','meg16_0327/161117','meg16_0330/161121','meg16_0332/161122','meg16_0333/161124','meg16_0337/161128','meg16_0339/161129','meg16_0340/161129','meg16_0341/161201','meg16_0343/161202','meg16_0345/161206','meg16_0346/161206','meg16_0348/161208','meg16_0349/161208','meg16_0350/161212','meg16_0352/161213'};
subjnum = [1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,6,2]; % counterbalancing numbers

%% parallelize
nw=32;
scheduler=cbu_scheduler();
scheduler.SubmitArguments='-q compute -l mem=170gb -l walltime=460800';
if isempty(gcp('nocreate')) || ~exist('pool','var') || pool.NumWorkers ~= nw,
    if ~isempty(gcp('nocreate'))
        delete(gcp('nocreate'))
    end
    scheduler.NumWorkers=nw;
    pool=parpool(scheduler,nw);
end

for sub = subs
    
    % go to subject directory
    cd(workingdir)
    swd = sprintf('sub%02d/%s',sub,subjects_dirs{sub}); % subject working directory
    cd(swd)
    
    %% load ROIs
    roi_folder = '/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp/ROIs';
%     roi_list = {'MD_ESV.nii','MD_IPS.nii','MD_AI.nii','MD_postMFG.nii','MD_ACCpreSMA.nii',...
%         'MD_FEF.nii','MD_antMFG.nii','MD_midMFG.nii','Auditory_Te3.nii','Visual_hOc1.nii'};
roi_list = {'LPFC.nii'};    
roinames=roi_list;
    roinames = strrep(roinames, '_', '-');
    roinames = strrep(roinames, '.nii', '');
    
    for rois = 1:numel(roi_list)
        Vroi = spm_vol(fullfile(roi_folder,roi_list{rois}));
        [Y, XYZ]=spm_read_vols(Vroi);
        coords = XYZ(:,Y>0)';
        
        %% attentional template
        % load data
        D = spm_eeg_load('caefMattn2_attention_task_block1_raw.mat');
        val = 1; %inversion method = MMN
        
        % mark button press trials
        bp_trials = []; % button press trials
        devents = D.events;
        for trial = 1:numel(D.events)
            if ismember(4096,[devents{trial}.value]) == 1
                bp_trials = [bp_trials trial];
            end
        end
        
        D = badtrials(D, bp_trials, 1);
        
        D.val=val;
        D.inv{val}.source.XYZ = coords;
        D.inv{val}.source.rad = 0;
        D.inv{val}.source.label= {roinames{rois}};
        D.inv{val}.source.fname = strcat(roinames{rois},'_attn2');
        D.inv{val}.source.type  = 'trials';
        spm_eeg_inv_extract_tw(D);
        
        Dattn = spm_eeg_load(strcat(roinames{rois},'_attn2'));
        data3D = single(Dattn.fttimelock.trial);
        %  sortedTNiNc = Dattn.fttimelock.trialinfo;
        conds = D.conditions;
        
        
        %% define labels
        
        % define subject specific pairings
        allperms = perms([1 2 3]);
        for i = 1:length(allperms)
            if rem(subjnum(sub),6)+1 == i
                choosetargets = allperms(i,:);
            end
        end
        target1 = sprintf('image%d',choosetargets(1));
        target2 = sprintf('image%d',choosetargets(2));
        target3 = sprintf('image%d',choosetargets(3));
        
        TNiNc = zeros(1,length(conds));
        indexcue1 = strfind(conds,'cue1');
        indexcue1 = ~cellfun(@isempty,indexcue1);
        for findcue1 = 1:length(indexcue1)
            if indexcue1(findcue1) == 1
                for i = 1:3 %three events per trial
                    try % use try because if the last cue does not follow three trials it will crash
                        if strcmp(conds(findcue1+i),'cue1')==1 || strcmp(conds(findcue1+i),'cue2')==1
                            break
                        elseif strcmp(conds(findcue1+i),target1) == 1
                            TNiNc(findcue1+i) = 0; %target (1item-cue1)
                        elseif strcmp(conds(findcue1+i),target2) == 1
                            TNiNc(findcue1+i) = 0; %inconsistent non-target (1item-cue1)
                        elseif strcmp(conds(findcue1+i),target3) == 1
                            TNiNc(findcue1+i) = 1; %consistent non-target (1item-cue1)
                        elseif regexp(conds{findcue1+i},'imageall') == 1
                            TNiNc(findcue1+i) = -1; %three-item display (3item-cue1)
                        else TNiNc(findcue1+i) = 0; %opaque
                        end
                    catch
                    end
                end
            end
        end
        indexcue2 = strfind(conds,'cue2');
        indexcue2 = ~cellfun(@isempty,indexcue2);
        for findcue2 = 1:length(indexcue2)
            if indexcue2(findcue2) == 1
                for i = 1:3 %three events per trial
                    try
                        if strcmp(conds(findcue2+i),'cue1')==1 || strcmp(conds(findcue2+i),'cue2')==1
                            break
                        elseif strcmp(conds(findcue2+i),target2) == 1
                            TNiNc(findcue2+i) = 0; %target (1item-cue2)
                        elseif strcmp(conds(findcue2+i),target1) == 1
                            TNiNc(findcue2+i) = 0; %inconsistent non-target (1item-cue2)
                        elseif strcmp(conds(findcue2+i),target3) == 1
                            TNiNc(findcue2+i) = 2; %consistent non-target (1item-cue2)
                        elseif regexp(conds{findcue2+i},'imageall') == 1
                            TNiNc(findcue2+i) = -2; %three-item display (3item-cue2)
                        else TNiNc(findcue2+i) = 0; %opaque
                        end
                    catch
                    end
                end
            end
        end
        
        % reformat to match Ds
        clabel = {};
        sortedTNiNc = [];
        trialtypes = D.inv{1}.inverse.trials;
        for i = 1:numel(trialtypes)
            keep = D.indtrial(trialtypes{i}, 'GOOD');
            clabel = [clabel D.conditions(keep)];
            sortedTNiNc = [sortedTNiNc TNiNc(keep)];
        end
        
        % define the two types of trials (with cue1 and cue2)
        type = cell(1,2);
        for i = 1:2
            type{i} = find(sortedTNiNc==i);
        end
        
        
        %% smoothing option
        smoothing = {'Raw','Smooth'};
        for smooth = 1:numel(smoothing)
            
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
                testwins      = [twin_test(1):stpsize:(twin_test(2)-winsize)]';
                testwins(:,2) = testwins(:,1) + winsize;
                ntests         = size(testwins,1);
                
                accuracy_matrix = nan(nfolds,ntrains,ntests);
                predicted_vals = cell(nfolds,ntrains,ntests);
                
                for fold = 1: nfolds
                    testind=[type{1}(fold:nfolds:end) type{2}(fold:nfolds:end)];
                    labels_test=sortedTNiNc(testind)';
                    
                    trainind=setdiff([type{:}],testind);
                    labels_train=sortedTNiNc(trainind)';
                    
                    parfor islide = 1:ntrains
                        itrain = data3D(:,:,trainwins(islide,1):trainwins(islide,2));
                        patternsTrain=double(itrain(trainind,:,:));
                        if smooth == 1 % no smoothing
                            patternsTrain=reshape(patternsTrain,length(trainind),[]);
                        elseif smooth == 2 % smoothing
                            patternsTrain = mean(patternsTrain,3);
                        end
                        
                        % classification using linear SVM (train the classifier)
                        model=svmtrain(labels_train,patternsTrain,libSVMsettings);
                        
                        for jslide = 1:ntests
                            itest = data3D(:,:,testwins(jslide,1):testwins(jslide,2));
                            patternsTest = double(itest(testind,:,:));
                            if smooth == 1 % no smoothing
                                patternsTest=reshape(patternsTest,length(testind),[]);
                            elseif smooth == 2 % smoothing
                                patternsTest = mean(patternsTest,3);
                            end
                            
                            % classification using linear SVM (test the classifier)
                            [predicted,accuracy,~]=svmpredict(labels_test,patternsTest,model);
                            accuracy_matrix(fold,islide,jslide)=accuracy(1);
                            predicted_vals{fold,islide,jslide} = predicted;
                            
                        end % next slide
                        
                    end % next slide
                end % next fold
                
                output_dir = 'crossgen_ROIs_1item_Nc_trialtarget';
                try cd(output_dir)
                catch eval(sprintf('!mkdir %s',output_dir)); cd(output_dir);
                end
                
                figure(sub)
                imagesc(squeeze(mean(accuracy_matrix,1)));
                colorbar;
                caxis([1,100]);
                
                % djm:
                set(gca,'xticklabel',(cellfun(@str2double,get(gca,'xticklabel'))-1)*stpsize*4-100);
                set(gca,'yticklabel',(cellfun(@str2double,get(gca,'yticklabel'))-1)*stpsize*4-100);
                xlabel('start of test windows (ms)')
                ylabel('start of train windows (ms)')
                
                saveas(gcf,sprintf('CrossGenMatrix%03dms%s_%s_sub%s.png',winsize*4,smoothing{smooth},roinames{rois},num2str(sub)));
                save(sprintf('CrossGenMatrix%03dms%s_%s_sub%s.mat',winsize*4,smoothing{smooth},roinames{rois},num2str(sub)));
                close;
                cd(workingdir)
                cd(swd)
                
            end % winsize
        end % smooth
    end % rois
    
end % sub

delete(gcp)