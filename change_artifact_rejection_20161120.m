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
subjects_dirs = {'meg16_0317/161107','meg16_0319/161110','meg16_0321/161111','meg16_0322/161114','meg16_0325/161115','meg16_0327/161117','meg16_0330/161121','meg16_0332/161122','meg16_0333/161124','meg16_0337/161128','meg16_0339/161129','meg16_0340/161129','meg16_0341/161201','meg16_0343/161202','meg16_0345/161206','meg16_0346/161206','meg16_0348/161208','meg16_0349/161208''meg16_0350/161212','meg16_0352/161213'};
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

for sub = [1,2,3,4,5,6,7,8,9,10,11,13]
    
    cd(out_pth)
    
    swd = sprintf('sub%02d/%s',sub,subjects_dirs{sub}); % subject working directory
    try cd(swd)
    catch eval(sprintf('!mkdir %s',swd)); cd(swd);
    end
    
    
    for run=1:subj_runs(sub) %length(block_names)
        
        %% Convert data
        
        for phase = 1:nphases(run)
            S = [];
            
            if phase == 1
                S = sprintf('efMspm12_%s_raw',block_names{run});
            elseif phase == 2
                S = sprintf('efMattn2_%s_raw',block_names{run});
            end

            
            D = spm_eeg_load(S);
            

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




