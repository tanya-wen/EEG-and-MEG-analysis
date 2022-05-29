%% Preprocessing of EEG + MEG data
% written by Tanya Wen
% 2016/09/21

% Modified from Script for MEEG preprocessing from Multimodal Example (SPM batch interface not yet sufficient) Rik Henson, July, 2012
%% Preprocessing of Preparatory Attention Study

addpath(genpath('/imaging/tw05'))
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'))
addpath(genpath('/imaging/local/software/spm_toolbox/eeglab13_4_3b'))
spm('defaults', 'eeg');

out_pth = '/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp'; % output path of processed data
raw_pth = '/megdata/cbu/prep_attention/';  % <--- INSERT PATH OF RAW DATA
cd(out_pth)

% Define BLOCK INFORMATION
block_names = {'auditory_template','visual_template','attention_task_block1','attention_task_block2','attention_task_block3','attention_task_block4','attention_task_block5'};
nphases = [1,1,2,2,2,2,2]; % attention task is cut into two phases for analysis
% Define EVENT INFORMATION
con_values{1} = [11,12,71,72]; % Trigger codes of interest (integers)
con_values{2} = [21,22,23,81,82,83];
con_values{3} = [1,2,31,32,33,34,35,36,37,38,39,91,92,93,94,95,96,97,98,99];
con_labels{1} =  {'cue1','cue2','cue1_mutate','cue2_mutate'}; % Condition labels corresponding to trigger codes
con_labels{2} =  {'image1','image2','image3','image1_opaque','image2_opaque','image3_opaque'};
con_labels{3} =  {'cue1','cue2','image1','image2','image3','imageall_pos1','imageall_pos2','imageall_pos3', ...
    'imageall_pos4','imageall_pos5','imageall_pos6','image1_opaque','image2_opaque','image3_opaque','imageall_pos1_opaque',...
    'imageall_pos2_opaque','imageall_pos3_opaque','imageall_pos4_opaque','imageall_pos5_opaque','imageall_pos6_opaque', };

% Define SUBJECT INFORMATION
subs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];  % subject numbers
subjects_dirs = {'meg16_0317/161107','meg16_0319/161110','meg16_0321/161111','meg16_0322/161114','meg16_0325/161115','meg16_0327/161117','meg16_0330/161121','meg16_0332/161122','meg16_0333/161124','meg16_0337/161128','meg16_0339/161129','meg16_0340/161129','meg16_0341/161201','meg16_0343/161202','meg16_0345/161206','meg16_0346/161206','meg16_0348/161208','meg16_0349/161208','meg16_0350/161212','meg16_0352/161213'};
% number of runs
subj_runs = 2+[4,5,5,5,5,5,5,4,5,5,5,5,5,5,5,5,5,5,5,5];
% User-defined bad EEG channels
bad_EEG = cell(1,max(subs));
bad_EEG{2}  = {'EEG021','EEG025'};
bad_EEG{3}  = {'EEG033'};
bad_EEG{4}  = {'EEG019','EEG073'};
bad_EEG{5}  = {'EEG002','EEG019','EEG070'};
bad_EEG{7}  = {'EEG002','EEG008','EEG065','EEG073'};
bad_EEG{8}  = {'EEG015','EEG021','EEG026','EEG073'};
bad_EEG{9}  = {'EEG018','EEG028','EEG029','EEG065','EEG008','EEG009','EEG016','EEG005','EEG035','EEG073'};
bad_EEG{10}  = {'EEG008','EEG028','EEG033','EEG038','EEG039','EEG054','EEG073','EEG065'};
bad_EEG{11}  = {'EEG039','EEG065','EEG073','EEG071','EEG074'};
bad_EEG{12}  = {'EEG027'};
bad_EEG{13}  = {'EEG029'};
bad_EEG{18} = {'EEG031'};
bad_EEG{19} = {'EEG039'};
bad_EEG{20} = {'EEG006','EEG029'};

chan = load(fullfile(out_pth,'meeg_channel_selection.mat'));

for sub = [1,2,3,4,5,6,7,8,9,10,11,13,15,16,17,18,19,20]
    
    cd(out_pth)
    
    swd = sprintf('sub%02d/%s',sub,subjects_dirs{sub}); % subject working directory
    try cd(swd)
    catch eval(sprintf('!mkdir %s',swd)); cd(swd);
    end
    
    
    for run=1:subj_runs(sub) %length(block_names)
        
        %% Convert data
        
        for phase = 1:nphases(run)
            S = [];
            S.dataset = fullfile(out_pth,swd,sprintf('%s_raw_ds.fif',block_names{run}));
            if phase == 1
                S.outfile = sprintf('spm12_%s_raw',block_names{run});
            elseif phase == 2
                S.outfile = sprintf('attn2_%s_raw',block_names{run});
            end
            S.channels = chan.label;
            S.save          = 0;
            S.reviewtrials  = 0;
            S.continuous    = 1;
            S.checkboundary = 0;
            
            D = spm_eeg_convert(S);
            
            
            %% User-defined bad EEG channels
            
            
            % Update user-defined bad EEG channels (not marked as bad in FIF file unfortunately)
            if ~isempty(bad_EEG{sub})
                D = badchannels(D,D.indchannel(bad_EEG{sub}),1);
            end
            
            % Downscale trigger channel so easier to see EOG when display "Other" channels
            ch = D.indchannel('STI101');
            trig_chan = D(ch,:,:);
            D(ch,:,:) = trig_chan/10^7;
            D.save;
            
            %% Could downsample at this point using spm_eeg_downsample, if space an issue
            
            
            
            %% ICA
            S = [];
            S.D = D;
            
            
            PCA_dim    = 60;   % number of PCs for data reduction prior to ICA (eg 60)
            VarThr     = 1/PCA_dim; % percentage of variance required to be an important artifactual IC (eg 1/PCA_dim, or could be 0) (DEFAULT possible)
            FiltPars   = [0.1,25,D.fsample];   % 1x2 or 1x3 matrix for low- or band- pass filtering before calculating correlation ([]=none=default), where third element is Sampling Rate
            TemAbsPval = 0.05; % p-value threshold for absolute temporal correlation to define artifact (eg 0.01)
            TemRelZval = 3;    % additional Z-value for relative threshold of absolute temporal correlation to define artifact (eg 3)
            SpaAbsPval = 0.05; % p-value threshold for absolute spatial correlation to define artifact (eg 0.01)
            SpaRelZval = 2;    % additional Z-value for relative threshold of absolute spatial correlation to define artifact (eg 2)
            PermPval   = 0.05; % p-value for boot-strapped absolute temporal correlation to define artifact (e.g, 0.05)
            Nperm      = 1000; % number of permutations to create null distribution of Pearson (eg 1000)
            
            modalities = {'MEGMAG';'MEGPLANAR';'EEG'};
            ref_chans = {'EOG061';'EOG062';'ECG063'};
            ref_type = {'VEOG','HEOG','ECG'};
            
            
            ica = {};
            for n = 1:numel(modalities)
                chans = indchantype(D,modalities{n},'GOOD'); % use data from good channels only for ica
                if PCA_dim > length(chans)
                    actual_PCA_dim = round(0.75*length(chans));
                    warning('Warning: For modality %s, ICA PCA dimensions set to %d (because only %d sensors)',modalitiles{n},actual_PCA_dim,length(chans));
                else
                    actual_PCA_dim = PCA_dim;
                end
                [weights,sphere,compvars,bias,signs,lrates,ICs] = runica(D(chans,:),'pca',actual_PCA_dim,'extended',1,'maxsteps',800); % will give different answer each time run
                ica{n}.weights = weights;
                ica{n}.chans = chans;
                ica{n}.compvars = compvars;
                ica{n}.modality = modalities{n};
            end
            clear weights chans % reuse variable names later in script
            
            % restrict to between first and last events to avoid modelling movemnts and other noise that we will not epoch
            triggevents = D.events;
            samp = round(triggevents(1).time*D.fsample) : round(triggevents(end).time*D.fsample);
            if rem(length(samp),2)==0  % Must be odd number of samples any fft below
                samp(1) = [];
            end
            Nsamp = size(samp,2);
            if rem(Nsamp,2)==0  % Must be odd number of samples for fft below
                Nsamp = Nsamp-1;
            end
            
            
            % Temporal reference signal for correlating with ICs
            for r = 1:length(ref_chans),
                refs.tem{r} = [];
                refs.tem{r} = D(indchannel(D,ref_chans{r}),samp);
            end
            if isempty(cell2mat(refs.tem)), warning('WARNING: No temporal reference has been found! Temporal artifacts cannot be removed!'); end
            
            % Spatial reference signal for correlating with ICs
            arttopos = load('/imaging/rh01/Methods/MEGEEGArtifactTemplateTopographies.mat');
            
            
            % Main loop
            PCA_dim = size(ica{1}.weights,1);
            toremove = cell(1,numel(ica)); TraMat  = cell(1,numel(ica));
            tempcor  = cell(1,numel(ica)); spatcor = cell(1,numel(ica));
            refs.spa = {}; artICs = {};
            for m = 1:numel(ica)
                
                modality{m}  = ica{m}.modality;
                weights{m}   = ica{m}.weights;
                chans{m}     = ica{m}.chans;
                
                iweights  = pinv(weights{m}); % weights will be in output file from previous step
                
                ICs = weights{m} * D(chans{m},samp);
                
                n = find(strcmp(arttopos.modalities,modality{m}));
                for r = 1:numel(ref_type)
                    if ~isfield(arttopos,ref_type{r})
                        warning('WARNING: No topograph found for %s, Spatial artifacts cannot be removed!',ref_type{r});
                        refs.spa{m}{r} = [];
                    else
                        tmp = getfield(arttopos,ref_type{r});
                        goodchans = ismember(indchantype(D,modalities{m}),indchantype(D,modalities{m},'GOOD'));
                        refs.spa{m}{r} = tmp{n}(goodchans)'; % detect_ICA_artefacts below assumes 1xChan vector
                    end
                end
                
                
                % Filtering (if any) (and transposition for speed)
                if length(FiltPars) == 3
                    fprintf('Bandpass filtering from %d to %d Hz (warning - can fail)',FiltPars(1), FiltPars(2));
                    for r=1:length(refs.tem)
                        refs.tem{r} = ft_preproc_bandpassfilter(refs.tem{r}, FiltPars(3), FiltPars(1:2),  [], 'but','twopass','reduce');
                    end
                    ICs  = ft_preproc_bandpassfilter(ICs,  FiltPars(3), FiltPars(1:2),  [], 'but','twopass','reduce')';
                elseif length(FiltPars) == 2
                    fprintf('Lowpass filtering to %d Hz',FiltPars(1));
                    for r=1:length(refs.tem)
                        refs.tem{r} = ft_preproc_lowpassfilter(refs.tem{r}, FiltPars(2), FiltPars(1),  5, 'but','twopass','reduce');
                    end
                    ICs  = ft_preproc_lowpassfilter(ICs,  FiltPars(2), FiltPars(1),  5, 'but','twopass','reduce')';
                else
                    ICs  = ICs'; % Faster for FFT  if pre-transpose once (check tic;toc)
                end
                
                
                % Check temporal correlation with any reference channels
                tempcor{m} = zeros(length(refs.tem),PCA_dim); tempval = zeros(length(refs.tem),PCA_dim);
                temprem  = cell(4,length(refs.tem));
                for r = 1:length(refs.tem)
                    if ~isempty(refs.tem{r})
                        
                        for k = 1:PCA_dim
                            [tempcor{m}(r,k),tempval(r,k)] = corr(refs.tem{r}',ICs(:,k));
                        end
                        
                        [~,temprem{1,r}] = max(abs(tempcor{m}(r,:)));
                        
                        temprem{2,r} = find(tempval(r,:) < TemAbsPval);
                        
                        temprem{3,r} = find(abs(zscore(tempcor{m}(r,:))) > TemRelZval);
                        
                        if Nperm > 0
                            permcor = zeros(1,PCA_dim);
                            maxcor  = zeros(Nperm,1);
                            
                            ff = fft(refs.tem{r}',Nsamp);
                            mf = abs(ff);
                            wf = angle(ff);
                            hf = floor((length(ff)-1)/2);
                            rf = mf;
                            
                            for l = 1:Nperm % could parfor...
                                rf(2:hf+1)=mf(2:hf+1).*exp((0+1i)*wf(randperm(hf)));    % randomising phases (preserve mean, ie rf(1))
                                rf((hf+2):length(ff))=conj(rf((hf+1):-1:2));            % taking complex conjugate
                                btdata = ifft(rf,Nsamp);                                % Inverse Fourier transform
                                
                                for k = 1:PCA_dim
                                    permcor(k) = corr(btdata,ICs(:,k));
                                end
                                maxcor(l) = max(abs(permcor));
                                fprintf('.');
                            end
                            fprintf('\n')
                            %         figure,hist(maxcor)
                            
                            temprem{4,r} = find(abs(tempcor{m}(r,:)) > prctile(maxcor,100*(1-PermPval)));
                        end
                    else
                        temprem{1,r} = []; temprem{2,r} = []; temprem{3,r} = []; temprem{4,r} = [];
                    end
                end
                
                
                % Check spatial correlation with any reference channels
                spatcor{m} = zeros(length(refs.spa{m}),PCA_dim); spatval = zeros(length(refs.spa{m}),PCA_dim);
                spatrem  = cell(3,length(refs.spa{m}));
                for r = 1:length(refs.spa{m})
                    if ~isempty(refs.spa{m}{r})
                        
                        for k = 1:PCA_dim
                            [spatcor{m}(r,k),spatval(r,k)] = corr(refs.spa{m}{r}',iweights(:,k));
                        end
                        
                        [~,spatrem{1,r}] = max(abs(spatcor{m}(r,:)));
                        
                        spatrem{2,r} = find(spatval(r,:) < SpaAbsPval);
                        
                        spatrem{3,r} = find(abs(zscore(spatcor{m}(r,:))) > SpaRelZval);
                    else
                        spatrem{1,r} = []; spatrem{2,r} = []; spatrem{3,r} = []; spatrem{4,r} = [];
                    end
                end
                
                
                % Variance Thresholding
                if VarThr > 0
                    compvars  = ica{m}.compvars;
                    varexpl   = 100*compvars/sum(compvars);
                    varenough = find(varexpl > VarThr);
                else
                    varenough = 1:PCA_dim;
                end
                
                
                % Define toremove
                for r = 1:length(refs.tem)
                    remove{m} = [1:PCA_dim];
                    for t = 1:3
                        remove{m} = intersect(remove{m},temprem{t,r});
                    end
                    if Nperm > 0; remove{m} = intersect(remove{m},temprem{4,r}); end
                    remove{m} = intersect(remove{m},varenough);  % plus sufficient Variance Explained (if required)
                    allrem{1,r} = remove{m};
                end
                
                for r = 1:length(refs.spa)
                    remove{m} = [1:PCA_dim];
                    for t = 1:3
                        remove{m} = intersect(remove{m},spatrem{t,r});
                    end
                    remove{m} = intersect(remove{m},varenough);  % plus sufficient Variance Explained (if required)
                    allrem{2,r} = remove{m};
                end
                
                all_remove{m} = [1:PCA_dim];
                if length(refs.tem{m}) > 0
                    all_remove{m} = intersect(all_remove{m},[allrem{1,:}]);
                end
                if length(refs.spa{m}) > 0
                    all_remove{m} = intersect(all_remove{m},[allrem{2,:}]);
                end
                
                if ~isempty(all_remove{m})
                    finalics  = setdiff([1:PCA_dim],all_remove{m});
                    TraMat{m} = iweights(:,finalics) * weights{m}(finalics,:);
                else
                    TraMat{m}    = eye(size(chans{m},2));
                end
                if m == 3; % for EEG
                    % Create the average reference montage
                    TraMat{m} = detrend(TraMat{m}, 'constant');
                end
            end
            
            
            %% Montage
            
            S = []; S.D = D; achans = [];
            for m = 1:numel(ica)
                achans = [achans chans{m}];
            end
            for c = 1:length(achans)
                S.montage.labelorg{c} = D.chanlabels{achans(c)};
            end
            S.montage.labelnew = S.montage.labelorg;
            S.montage.tra      = blkdiag(TraMat{:});
            S.keepothers       = 1;
            S.prefix           = 'M';
            D = spm_eeg_montage(S);
            
            %% Bandpass filter
            S = [];
            S.D = D.fname;
            S.type = 'butterworth';
            S.order = 5;
            S.band = 'bandpass';
            S.freq = [0.1 40];
            D = spm_eeg_filter(S);
            
            
            %% Define epochs
            S = [];
            S.D = D;
            
            if run == 1 % auditory template
                S.timewin = [-100 1000];
                for n = 1:length(con_labels{run})
                    S.trialdef(n).conditionlabel = con_labels{run}{n};
                    S.trialdef(n).eventtype = 'STI101_up';
                    S.trialdef(n).eventvalue = con_values{run}(n);
                    S.trialdef(n).trlshift = 26;
                end
            end
            if run == 2 % visual template
                S.timewin = [-100 1500]; 
                for n = 1:length(con_labels{run})
                    S.trialdef(n).conditionlabel = con_labels{run}{n};
                    S.trialdef(n).eventtype = 'STI101_up';
                    S.trialdef(n).eventvalue = con_values{run}(n);
                    S.trialdef(n).trlshift = 34;
                end
            end
            if run >= 3 % attention tasks
                if phase == 1 % preparatory attention phase
                    S.timewin = [-100 1750];
                    for n = 1:length(con_labels{3})
                        S.trialdef(n).conditionlabel = con_labels{3}{n};
                        S.trialdef(n).eventtype = 'STI101_up';
                        S.trialdef(n).eventvalue = con_values{3}(n);
                        S.trialdef(n).trlshift = 26;
                    end
                elseif phase == 2 % stimulus encoding phase
                    S.timewin = [-100 1500];
                    for n = 1:length(con_labels{3})
                        S.trialdef(n).conditionlabel = con_labels{3}{n};
                        S.trialdef(n).eventtype = 'STI101_up';
                        S.trialdef(n).eventvalue = con_values{3}(n);
                        S.trialdef(n).trlshift = 34;
                    end
                end
            end
            
            S.bc = 1;
            S.prefix = 'e';
            S.eventpadding = 0;
            D = spm_eeg_epochs(S);
            
            
            %% Artifact
            S = [];
            S.D = D;
            S.methods(1).fun = 'flat';
            S.methods(1).channels = 'MEG';
            S.methods(1).settings.threshold = 0;
            S.methods(1).settings.seqlength = 2;
            S.methods(end+1).fun = 'flat';
            S.methods(end).channels = 'EEG';
            S.methods(end).settings.threshold = 0;
            S.methods(end).settings.seqlength = 2;
            S.methods(end+1).fun = 'peak2peak';
            S.methods(end).channels = 'MEGMAG';
            S.methods(end).settings.threshold = 4000;
            S.methods(end+1).fun = 'peak2peak';
            S.methods(end).channels = 'MEGPLANAR';
            S.methods(end).settings.threshold = 400;
            S.methods(end+1).fun = 'peak2peak';
            S.methods(end).channels = 'EEG';
            S.methods(end).settings.threshold = 120;
            S.methods(end+1).fun = 'peak2peak';
            S.methods(end).channels = 'EOG';
            S.methods(end).settings.threshold = 750;
            
            D = spm_eeg_artefact(S);
            
            
            %% Remove response trials
            
            % mark button press trials
            bp_trials = []; % button press trials
            devents = D.events;
            for trial = 1:numel(D.events)
                if ismember(4096,[devents{trial}.value]) == 1
                    bp_trials = [bp_trials trial];
                end
            end
      
            D = badtrials(D, bp_trials, 1);
 
            
            %% Averaging
            
            S   = [];
            S.D = D;
            S.robust.ks = 3;
            S.robust.bycondition = false;
            S.robust.savew = false;
            S.robust.removebad = true;
            S.circularise = false;
            S.prefix = 'm';
            D = spm_eeg_average(S);
            
        end % number of phases
        
    end % number of runs
    
end % number of subjects

cd(out_pth)




