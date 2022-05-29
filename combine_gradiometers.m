%% Combined gradiometers

clear all
clc

addpath /usr/local/spm12
addpath /work/imaging3/MEG/Camcan_word_recognition/scripts/spm_melek

spm('defaults', 'EEG')

%% Define your subject data MEG
% subjects option

data_path = 
'/work/imaging3/MEG/Camcan_word_recognition/data/data_B_trans/aamod_meg_denoise_ICA_2_applytrajectory_00002/';
addpath 
/work/imaging3/MEG/Camcan_word_recognition/functions/QueryFunction/QueryFun_v1;

q = [];
q.SessionList = {
     'MEG' 
'/work/imaging3/MEG/Camcan_word_recognition/data/data_B_trans/aamod_meg_denoise_ICA_2_applytrajectory_00002/*<ccID>/analysis_orthography_2/bdef_Mspm12_transdef_transrest_mf2pt2_wordrecog_raw.mat';
     };
q = CCQuery_CheckFiles(q);

info = LoadSubIDs;
CCIDS = info.SubCCIDc;

%% options

S = [];
S.mode = 'scalp x time';
S.conditions = {'consonants';'words'};
S.timewin = [-100 600];

% dobascor = 1;


%% Loop through subjects
for idx = 1:numel(q.AllExistSubInd)
     ii = q.AllExistSubInd(idx);

     fname = q.FileNames.MEG{ii};
     subjname = splitstring(fname,'/');
     subjname = subjname{8};

     this_dir = [data_path subjname '/' 'analysis_orthography_2'];


%     this_megcomb = [this_dir 
'/PPpmbdef_Mspm12_trans_mf2pt2_wordrecog_raw.mat'];
%     if exist (this_megcomb, 'file') == 0 % if processing didn't run 
yet

         fprintf('\n**** Processing subject %s ****\n\n',CCIDS{ii})

         % write file

         input_MEG = [this_dir 
'/bdef_Mspm12_transdef_transrest_mf2pt2_wordrecog_raw.mat'];

         SS.D = input_MEG;
         SS.mode = 'replacemeg';
         SS.prefix = 'P';
         D = spm_eeg_combineplanar(SS);

%         clear D;
%         S = [];
%         S.D = [this_dir 
'/Pbdef_Mspm12_transdef_transrest_mf2pt2_wordrecog_raw.mat'];
%         S.timewin = [-100 0];
%         S.prefix='b';
%         D = spm_eeg_bc(S);

end
