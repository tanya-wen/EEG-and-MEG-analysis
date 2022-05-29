%% Cross-generalization of auditory localizer to preparatory attention task. Done separately within each ROI
% written by Tanya Wen
% 2017/01/26

dbstop if error

addpath('/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp')
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
%addpath(genpath('/imaging/local/software/spm_toolbox/eeglab13_4_3b'))
spm('defaults', 'eeg');

% required software: libSVM & RSAtoolbox
workingdir = '/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp';
addpath(fullfile('/imaging/tw05/Preparatory_Attention_Study/','software','libsvm-mat-2.87-1'));
% control variables
libSVMsettings='-s 1 -t 0'; % nu-SVM, linear
nRandomisations=1000;
% rmpath('/hpc-software/matlab/r2014a/toolbox/bioinfo/biolearning/'); % to make sure libSVM code is used (not strictly necessary: matlab svmtrain yields exactly same model)

% Define SUBJECT INFORMATION
subs = [1,2,3,4,5,6,7,8,9,10,11,13,15,16,17,18,19,20];  % subject numbers
subjects_dirs = {'meg16_0317/161107','meg16_0319/161110','meg16_0321/161111','meg16_0322/161114','meg16_0325/161115','meg16_0327/161117','meg16_0330/161121','meg16_0332/161122','meg16_0333/161124','meg16_0337/161128','meg16_0339/161129','meg16_0340/161129','meg16_0341/161201','meg16_0343/161202','meg16_0345/161206','meg16_0346/161206','meg16_0348/161208','meg16_0349/161208','meg16_0350/161212','meg16_0352/161213'};
subjnum = [1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,6,2]; % counterbalancing numbers

%% parallelize
nw=120;
scheduler=cbu_scheduler();
scheduler.SubmitArguments='-q compute -l mem=300gb -l walltime=460800';

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
    

    %% visual template
    % load data
    D = spm_eeg_load('aefMspm12_visual_template_raw.mat');
    data3Dtest = D.fttimelock.trial;
    
    % mark button press trials
    bp_trials = []; % button press trials
    devents = D.events;
    for trial = 1:numel(D.events)
        if ismember(4096,[devents{trial}.value]) == 1
            bp_trials = [bp_trials trial];
        end
    end
    try
    D = badtrials(D, bp_trials, 1);
    catch
    end
    
    % get good channels
    modalities = {'MEGMAG';'MEGPLANAR';'EEG'};
    for m = 1:numel(modalities)
        chans1{m} = indchantype(D,modalities{m},'GOOD');
    end
    
    cond_values_test = D.fttimelock.trialinfo;
    
    % define labels
    
    % define image1, image2, and image3 index vectors
    images=cell(1,3);
    for i=1:3
        images{i} = D.indtrial(sprintf('image%d',i), 'GOOD');
    % images{i} = find(cond_values_train==i);
    end
    
    testind=[images{1} images{2} images{3}];
    labels_test=cond_values_test(testind)';
    
    clear D
    
    %% attentional template
    % load data
    D = spm_eeg_load('caefMspm12_attention_task_block1_raw.mat');
    data3Dtrain = D.fttimelock.trial;
    
    % mark button press trials
    bp_trials = []; % button press trials
    devents = D.events;
    for trial = 1:numel(D.events)
        if ismember(4096,[devents{trial}.value]) == 1
            bp_trials = [bp_trials trial];
        end
    end
    
    try
        D = badtrials(D, bp_trials, 1);
    catch
    end
    
    % get good channels
    modalities = {'MEGMAG';'MEGPLANAR';'EEG'};
    for m = 1:numel(modalities)
        chans2{m} = indchantype(D,modalities{m},'GOOD');
    end
    
    cond_values_train = D.fttimelock.trialinfo;

   
    % define labels

    % define subject specific pairings
    allperms = perms([1 2 3]);
    for i = 1:length(allperms)
        if rem(subjnum(sub),6)+1 == i
            choosetargets = allperms(i,:);
        end
    end
    
    % define cue1 & cue2 index vectors
    cues=cell(1,2);
    for i=1:2
        cues{i} = D.indtrial(sprintf('cue%d',i), 'GOOD');
        % cues{i} = find(cond_values_test==i);
    end
    
    trainind=[cues{1} cues{2}];
    labels_train=cond_values_train(trainind)';
    labels_train(labels_train==1)=-1;
    labels_train(labels_train==2)=-2;
    labels_train(labels_train==-1)=choosetargets(1);
    labels_train(labels_train==-2)=choosetargets(2);
    
    % djm:
    skiptest=~ismember(labels_test,unique(labels_train));
    testind(skiptest)=[];
    labels_test(skiptest)=[];
    
    newdata3Dtrain = zscore(data3Dtrain(:,intersect([chans1{:}],[chans2{:}]),:),[],2); % standardize each channel
    newdata3Dtest = zscore(data3Dtest(:,intersect([chans1{:}],[chans2{:}]),:),[],2); % standardize each channel
    
    clear D
    
    %% smoothing option
    smoothing = {'Raw','Smooth'};
    for smooth = 1:numel(smoothing)
        
        %% sliding time window
        winsizes = 8;%[8,25,50,125]; %32ms,100ms,200ms,500ms
        for winsize = winsizes
            
            % set analysis time windows
            twin_train = [1 size(newdata3Dtrain,3)]; % time of interest
            stpsize        = 1; %step size
            trainwins      = [twin_train(1):stpsize:(twin_train(2)-winsize)]';
            trainwins(:,2) = trainwins(:,1) + winsize;
            ntrains        = size(trainwins,1);
            twin_test = [1 size(newdata3Dtest,3)]; % time of interest
            testwins      = [twin_test(1):stpsize:(twin_test(2)-winsize)]';
            testwins(:,2) = testwins(:,1) + winsize;
            ntests         = size(testwins,1);
            
            accuracy_matrix = nan(ntrains,ntests);
            predicted_vals = cell(ntrains,ntests);
            parfor islide = 1:ntrains
                itrain = newdata3Dtrain(:,:,trainwins(islide,1):trainwins(islide,2));
                patternsTrain=double(itrain(trainind,:,:));
                if smooth == 1 % no smoothing
                    patternsTrain=reshape(patternsTrain,length(trainind),[]);
                elseif smooth == 2 % smoothing
                    patternsTrain = mean(patternsTrain,3);
                end
                % classification using linear SVM (train the classifier)
                model=svmtrain(labels_train,patternsTrain,libSVMsettings);
                
                for jslide = 1:ntests
                    itest = newdata3Dtest(:,:,testwins(jslide,1):testwins(jslide,2));
                    patternsTest = double(itest(testind,:,:));
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
            
            output_dir = 'crossgen_sensor_prep-vis';
            try cd(output_dir)
            catch eval(sprintf('!mkdir %s',output_dir)); cd(output_dir);
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
            
            saveas(gcf,sprintf('CrossGenMatrix%03dms%s_sub%s.png',winsize*4,smoothing{smooth},num2str(sub)));
            save(sprintf('CrossGenMatrix%03dms%s_sub%s.mat',winsize*4,smoothing{smooth},num2str(sub)));
            close;
            cd(workingdir)
            cd(swd)
            
        end % winsize
    end % smooth
    
end

delete(gcp)