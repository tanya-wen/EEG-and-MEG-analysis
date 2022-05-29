%% Plotting topography over a period of interest 
% Plots scalp data averaged over a period of time (e.g. significant cluster % from statistical testing

clear

addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
addpath('/imaging/local/software/spm_toolbox/eeglab13_4_3b')
spm('defaults', 'eeg');

workingdir = '/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp';

% Define SUBJECT INFORMATION
subs = [1,2,3,4,5,6,7,8,9,10,11,13,15,16,17,18,19,20];  % subject numbers
subjects_dirs = {'meg16_0317/161107','meg16_0319/161110','meg16_0321/161111','meg16_0322/161114','meg16_0325/161115','meg16_0327/161117','meg16_0330/161121','meg16_0332/161122','meg16_0333/161124','meg16_0337/161128','meg16_0339/161129','meg16_0340/161129','meg16_0341/161201','meg16_0343/161202','meg16_0345/161206','meg16_0346/161206','meg16_0348/161208','meg16_0349/161208','meg16_0350/161212','meg16_0352/161213'};
subjnum = [1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,6,2]; % counterbalancing numbers
subj_eeg = [1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1]; % whether to analyze EEG data


task = 'mcaefMspm12_attention_task_block1_raw.mat';
% task = 'mN2pc_caefMattn2_attention_task_block1_raw.mat';


%% calculate Grand Mean across subjects
fname = {};
for sub = subs
%     if subj_eeg(subs) ~= 0
    cd(workingdir)
    swd = sprintf('sub%02d/%s',sub,subjects_dirs{sub}); % subject working directory
    cd(swd)
    fname{end+1} = fullfile(workingdir,swd,task);
%     end
end

cd(workingdir)
S.D = char(fname);
S.outfile = sprintf('new_grand_mean_%s', task);
%D = spm_eeg_grandmean_tw(S);


%% plot topography from t1 to t2

t1 = 250; % start [ms]
t2 = 300; % end [ms]

D = spm_eeg_load(sprintf('/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp/sub01/meg16_0317/161107/new_grand_mean_%s', task)); % change for grad/mag

three_item_cond = find(~cellfun(@isempty,regexp(D.conditions,'imageall_pos\d$')));



%%%%%
try
    addpath /imaging/dm01/MEG/aaMEG/
    [LHS, RHS, trueloc, midline, MAGS, LONGS, LATS]=getLRpairs2(D);
catch
    x=load('/imaging/dm01/MEG/aaMEG/LRpairs.mat');
end


modalities = {'EEG','MAGS','LONGS','LATS'};            
for m = 1:numel(modalities)
    in.f=m+1;
    try delete(in.f);
    catch
    end
    figure(in.f); in.noButtons = true; set(in.f,'color','w');
    switch modalities{m}
        case 'EEG'
            chanind=D.indchantype('EEG');
        otherwise
            chanind=x.(modalities{m});
    end
    [qZ, qf, qd]=spm_eeg_plotScalpData_d(...
        squeeze(mean(mean(D(chanind,D.indsample(1e-3*t1):D.indsample(1e-3*t2), three_item_cond), 2),3)), ...
        D.coor2D(chanind),D.chanlabels(chanind),in);
    
    fpos=D.coor2D(chanind);
    xmin    = min(fpos(1,:));
    xmax    = max(fpos(1,:));
    dx      = (xmax-xmin)./100;
    ymin    = min(fpos(2,:));
    ymax    = max(fpos(2,:));
    dy      = (ymax-ymin)./100;

end
%%%%%%
