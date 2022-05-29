%% Source localization of EEG + MEG data
% written by Tanya Wen
% 2016/09/23

%% Source localization of Preparatory Attention Study

addpath(genpath('/imaging/tw05'))
% addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'))
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6685'))
addpath(genpath('/imaging/local/software/spm_toolbox/eeglab13_4_3b'))
spm('defaults', 'eeg');

out_pth = '/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp'; % output path of processed data
raw_pth = '/megdata/cbu/prep_attention/';  % <--- INSERT PATH OF RAW DATA
cd(out_pth)

% Define SUBJECT INFORMATION
subs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];  % subject numbers
subjects_dirs = {'meg16_0317/161107','meg16_0319/161110','meg16_0321/161111','meg16_0322/161114','meg16_0325/161115','meg16_0327/161117','meg16_0330/161121','meg16_0332/161122','meg16_0333/161124','meg16_0337/161128','meg16_0339/161129','meg16_0340/161129','meg16_0341/161201','meg16_0343/161202','meg16_0345/161206','meg16_0346/161206','meg16_0348/161208','meg16_0349/161208','meg16_0350/161212','meg16_0352/161213'};
% number of runs
subj_runs = [4,5,5,5,5,5,5,4,5,5,5,5,5,5,5,5,5,5,5,5];
subj_eeg = [1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1]; % whether to analyze EEG data

fnam = {'aefMspm12_auditory_template_raw.mat','aefMspm12_visual_template_raw.mat','caefMspm12_attention_task_block1_raw.mat','caefMattn2_attention_task_block1_raw.mat'};

display_flag = 1;

for sub = [9,10,11,13]
    
    cd(out_pth)
    swd = sprintf('sub%02d/%s',sub,subjects_dirs{sub}); % subject working directory
    
    cd(swd)
    
    for file=1:length(fnam)
        %% merge data
        if strcmp(fnam{file}, 'caefMspm12_attention_task_block1_raw.mat')==1
            S.D = [];
            for block = 1:subj_runs(sub) % the five attention blocks
                S.D = strvcat(S.D,sprintf('aefMspm12_attention_task_block%s_raw.mat',num2str(block)));
            end
            S.recode='same';
            D = spm_eeg_merge(S);
        end
        if strcmp(fnam{file}, 'caefMattn2_attention_task_block1_raw.mat')==1
            S.D = [];
            for block = 1:subj_runs(sub) % the five attention blocks
                S.D = strvcat(S.D,sprintf('aefMattn2_attention_task_block%s_raw.mat',num2str(block)));
            end
            S.recode='same';
            D = spm_eeg_merge(S);
        end
        
        %% Create Forward models
        S=[]; S.D = fullfile(out_pth,swd,fnam{file});
        D = spm_eeg_load(S.D);
        
        val = 1;
        D.inv{val} = [];  % If want to be safe!
        D.val = val;
        D.inv{val}.date    = strvcat(date,datestr(now,15));
        D.inv{val}.comment = 'Standard';
        
        %%% Mesh (if spm_eeg_inv_mesh were more flexible, then could pick up normalisation
        %%% params if calculated already (for fMRI analysis)!)
        T1nam = dir(fullfile(out_pth,sprintf('sub%02d/%s',sub),'MPRAGE','co*.nii'));
        D.inv{val}.mesh.sMRI = fullfile(out_pth,sprintf('sub%02d/%s',sub),'MPRAGE',T1nam.name);
        D = spm_eeg_inv_mesh_ui(D, val, D.inv{val}.mesh.sMRI, 2);
        
        %%% Data Reg
        fids = load(fullfile(out_pth,sprintf('sub%02d/%s',sub),'MPRAGE',sprintf('mri_fids_sub%02d.mat',sub)));
        MRIfids           = [];
        MRIfids.fid.pnt   = fids.mri_fids;
        MRIfids.fid.label = {'Nasion';'LPA';'RPA'};
        MRIfids.pnt       = D.inv{val}.mesh.fid.pnt;        % Scalp mesh points from MRI above
        MRIfids.unit      ='mm';
        MEGfids           = D.fiducials;  % To remove nose points in next line
        MEGfids.pnt       = MEGfids.pnt(find(~(MEGfids.pnt(:,2)>0 & MEGfids.pnt(:,3)<0)),:);
        
        D = spm_eeg_inv_datareg_ui(D, val, MEGfids, MRIfids, 1);
        
        %%% Forward modelling
        D.inv{val}.forward = struct([]);
        for ind = 1:length(D.inv{val}.datareg)
            if strcmp(D.inv{val}.datareg(ind).modality,'EEG')
                D.inv{val}.forward(ind).voltype = 'EEG BEM';
                %            D.inv{val}.forward(ind).voltype = 'OpenMEEG BEM';
            elseif strcmp(D.inv{val}.datareg(ind).modality,'MEG')
                D.inv{val}.forward(ind).voltype = 'Single Shell';
                %            D.inv{val}.forward(ind).voltype = 'Single Sphere';
            end
        end
        D = spm_eeg_inv_forward(D);
        
        D.save
        
        
        %% Display Results
        
        if display_flag
            
            cd(out_pth)
            swd = sprintf('sub%02d/%s',sub,subjects_dirs{sub}); % subject working directory
            
            cd(swd)
            S=[]; S.D = fullfile(out_pth,swd,fnam{file});
            D = spm_eeg_load(S.D);
            
            spm_eeg_inv_checkmeshes(D);
            
            for ind = 1:length(D.inv{1}.datareg)
                spm_eeg_inv_checkdatareg(D, 1, ind);
            end
            
            for ind = 1:length(D.inv{1}.datareg)
                spm_eeg_inv_checkforward(D, 1, ind);
            end
            
        end
        
        %% Do basic MMN and MSP inversions
        % MMN: minimum norm estimation; MSP: multiple sparse priors; EBB: Empirical Bayes Beamformer
        if file == 1 || file == 3
            inv_conditions = {'cue1','cue2'};  % Could localise difference
        elseif file == 2
            inv_conditions = {'image1','image2','image3'};
        elseif file == 4
            inv_conditions = {'image1','image2','image3','imageall_pos1','imageall_pos2','imageall_pos3','imageall_pos4','imageall_pos5','imageall_pos6'};
        end
        if  subj_eeg(sub) == 1 %add this part due to EEG problems
            inv_modalities = {{'MEGMAG' 'MEGPLANAR' 'EEG'} {'MEGMAG' 'MEGPLANAR' 'EEG'} {'MEGMAG' 'MEGPLANAR' 'EEG'}} % 'MEGMAG' 'MEGPLANAR' 'EEG' };
        elseif subj_eeg(sub) == 0 %add this part due to EEG problems
            inv_modalities = {{'MEGMAG' 'MEGPLANAR'} {'MEGMAG' 'MEGPLANAR'} {'MEGMAG' 'MEGPLANAR'}} % 'MEGMAG' 'MEGPLANAR' 'EEG' };
        end
        inv_type =       {'MMN', 'GS', 'EBB'} % 'MMN'    'MMN'    'MMN'};
        inv_twin = [-100 4000];   % Late familiarity effect may be removed by Hanning
        inv_fboi = [0 40];    % All since data already filtered
        
        F = []; Nm = []; VarExp = [];
        cd(out_pth)
        
        
        cd(out_pth)
        swd = sprintf('sub%02d/%s',sub,subjects_dirs{sub}); % subject working directory
        
        cd(swd)
        S=[]; S.D = fullfile(out_pth,swd,fnam{file});
        
        D.inv{val}.inverse = [];   % Clear to be safe!
        
        val = 0;
        for m = 1:length(inv_type);   % First without fMRI priors
            val = val+1;
            D.val = val;
            D.inv{val}.comment = {sprintf('Ind: Val%d: Mod:%s Inv:%s',val,cat(2,inv_modalities{m}{:}),inv_type{m})};
            D.inv{val}.mesh    = D.inv{1}.mesh;
            D.inv{val}.datareg = D.inv{1}.datareg;
            D.inv{val}.forward = D.inv{1}.forward;
            try D.inv{val}.gainmat = D.inv{1}.gainmat; end  % Doesn't exist on very first inversion
            
            D.inv{val}.inverse.trials   = inv_conditions;
            D.inv{val}.inverse.woi      = inv_twin;
            D.inv{val}.inverse.Han      = 0;  %For late Fam effect???
            D.inv{val}.inverse.lpf      = inv_fboi(1);
            D.inv{val}.inverse.hpf      = inv_fboi(2);
            D.inv{val}.inverse.type     = inv_type{m};
            D.inv{val}.inverse.modality = inv_modalities{m};
            D = spm_eeg_invert(D);
            
            F(sub,val)      = D.inv{val}.inverse.F;
            Nm(sub,val)     = size(D.inv{val}.inverse.L,1);
            VarExp(sub,val) = D.inv{val}.inverse.VE;
        end
        
        D.save;
        
        
        %% Review results
        if display_flag
            
            cd(out_pth)
            swd = sprintf('sub%02d/%s',sub,subjects_dirs{sub}); % subject working directory
            
            cd(swd)
            S=[]; S.D = fullfile(out_pth,swd,fnam{file});
            D = spm_eeg_load(S.D);
            
            for v=1:length(D.inv)
                D.val=v;
                spm_eeg_invert_display(D,170);
                disp([sub v])
            end
            
        end
    end
end

%% Write time-freq contrasts for all inversions above
Src_images={};
cvals = [1 2 3];
con_twin = inv_twin;
%con_twin = [140 200]; % for just N170?
con_fboi = inv_fboi;
con_engy = 'evoked';

for sub = subs
    
    swd = sprintf('sub%02d/%s',sub,subjects_dirs{sub}); % subject working directory
    
    for file=1:length(fnam)
        D = spm_eeg_load(fullfile(out_pth,swd,fnam{file}));
        % Choose best inversion (first here?)
        
        for val = cvals
            D.val = val;
            
            % Calculate a time-freq contrast...
            
            D.inv{val}.contrast.woi  = con_twin;
            D.inv{val}.contrast.fboi = con_fboi;
            D.inv{val}.contrast.type = con_engy;
            D = spm_eeg_inv_results(D);
            
            %...and write it as an image
            
            D.inv{val}.contrast.display   = 0;
            D.inv{val}.contrast.space     = 1;  % MNI
            D.inv{val}.contrast.smoothing = 12;  % Could play with (8=default)
            D = spm_eeg_inv_Mesh2Voxels(D);
            
            Src_images{val}{sub} = D.inv{val}.contrast.fname;
        end
        
    end
end