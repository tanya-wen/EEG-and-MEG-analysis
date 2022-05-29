  % Maxfilter 2.2 Matlab script
%
% script adapted from Olaf's Maxfilter_main_script_OH.m
% LJM January 2016
%
% This version has three options for running on our Compute machines, based
% on setting of ParType variable below: 0 = run on Login in serial; 1 = run
% one Compute core (using spmd); 2 = run on multiple Compute cores
% (using parfor) - with one core per subject 
%
% Maxfilter is applied in two steps (see commands below):
%   1. SSS with autobad on, to detect bad channels
%   2. tSSS, downsampling, and interpolate to first file if multiple runs
%
% MF command as part of bad channel detection:
% !/neuro/bin/util/x86_64-pc-linux-gnu/maxfilter-2.2  -f /data_in_dir/sub_dirs/raw_files_in.fif -o /data_out_dir/sub_dirs/raw_files_out_bad.fif 
% -origin x y z -frame head -ctc /neuro/databases/ctc/ct_sparse.fif -cal /neuro/databases/sss/sss_cal.dat -linefreq 50  -force -autobad 900 -badlimit 7  
% -headpos -hpistep 10 -hp /data_out_dir/sub_dirs/raw_files_out_headpos.txt -v | tee /data_out_dir/sub_dirs/raw_files_out_bad.log
%
% Final MF command:
% !/neuro/bin/util/x86_64-pc-linux-gnu/maxfilter-2.2  -f /data_in_dir/sub_dirs/raw_files_in.fif -o /data_out_dir/sub_dirs/raw_files_out_ssst.fif 
% -ctc /neuro/databases/ctc/ct_sparse.fif -cal /neuro/databases/sss/sss_cal.dat -linefreq 50  -force -autobad 900 -badlimit 7   -st 10 -corr 0.98 
% -movecomp inter  -hpisubt amp -origin x y z -frame head -trans /data_out_dir/sub_dirs/raw_files_out_ssst_trans1st.fif  -ds 4 -v 
% | tee /data_out_dir/sub_dirs/raw_files_out_ssst_trans1st.log
%
% MF 2.2 manual:
% http://imaging.mrc-cbu.cam.ac.uk/meg/Maxfilter_V2.2?action=AttachFile&do=view&target=Maxfilter_Manual_v2pt2.pdf
%
% Olaf Hauk, Sep 2014 (modified from RH's script)

addpath /imaging/local/meg_misc
addpath /neuro/meg_pd_1.2/

addpath( pwd ) % for local functions

%% General options
% 0: Run on Login machines (frowned upon!)
% 1: Run maxfilter call on Compute machines using spmd (faster)
% 2: Run on multiple Compute machines using parfor (best, but less feedback if crashes)
ParType = 0;   % parallelisation options 

% overwrite exisiting fiff-files?
datinfo.OverWrite = 1;  % Set to 1 if you want to overwrite existing files

% if you don't want to run maxfilter, but just look at command line output,
% then specify do_eval = 0
datinfo.do_eval = 1;

% which subjects from list below to process ([] if all)
datinfo.subj_list = [];
% which runs per subject to process ([] if all)
datinfo.run_list{1} = [];

% general log-file, records which fiff-files have been successfully created
% will be created for each subject, beginning with subject index
% separate files for different subjects because of parfor-loop
datinfo.logfilestem = 'general_log_tmp.log';

%% Data info - needs to be edited for every study
% Raw data will be read from subdirectories of...
datinfo.data_in_dir = '/megdata/cbu/prep_attention/'; 

% Results will be written to subdirectories of...
datinfo.data_out_dir = '/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp/sub20';

% sub-directories for individual subjects
datinfo.sub_dirs = {'meg16_0352/161213'};


% input names of raw data files - can be different for different subjects, if only one list provided, it will be used for all subjects
datinfo.raw_files_in = {{'auditory_template_raw', 'visual_template_raw', 'attention_task_block1_raw', 'attention_task_block2_raw', 'attention_task_block3_raw', 'attention_task_block4_raw','attention_task_block5_raw'}};

% output names of raw data files - if only one list provided, it will be used for all subjects
datinfo.raw_files_out = {{'auditory_template_raw', 'visual_template_raw', 'attention_task_block1_raw', 'attention_task_block2_raw', 'attention_task_block3_raw', 'attention_task_block4_raw','attention_task_block5_raw'}};

% specify files for which movement compensation shall NOT be applied
% put nothing "{}" if mvcomp SHALL be applied to all files 
datinfo.mvcomp_fail = {};

% bad channels, per subject in sequence of cbu_codes{}
% if no bad channels, add {''} for that subject
datinfo.bad_chans = {{'1943','0813'}};

datinfo.del_bads = 1;   % 1: delete fiff-files created during bad channel detection; 0: don't

%% general parameters (may need editing)

% info from movement compensation
datinfo.movfile = 'trans_move.txt'; % This file will record translations between runs

% Number of data buffers for each experiment
datinfo.Nbuf = [480];   


%% Maxfilter parameters
% Maxfilter command string options in structure format (will be concatenated)
% YOU MAY WANT TO EDIT SOME OF THIS
cstr.dsstr = ' -ds 4 ';   % DOWNSAMPLE from 1000 to 250Hz
%cstr.dsstr = '';
cstr.headori = [];  % if empty, sphere will be fitted to the digi points
%cstr.maxfstr = '!/neuro/bin/util/x86_64-pc-linux-gnu/maxfilter-2.2 '; % this specifies the Maxfilter version
cstr.maxfstr = '!/hpc-software/neuromag/bin/util/maxfilter-2.2.12 ';
%cstr.maxfstr = '!maxfilter-2.2.12 ';
basestr = ' -ctc /neuro/databases/ctc/ct_sparse.fif -cal /neuro/databases/sss/sss_cal.dat';
basestr = [basestr ' -linefreq 50  -force'];
cstr.basestr = basestr;
cstr.badstr  = sprintf(' -autobad %d -badlimit %d ',900,7); % 900s is 15mins - ie enough for whole recording
cstr.compstr_bad = ' -headpos -hpistep 10';  % for bad channel detection only
cstr.compstr = sprintf(' -movecomp inter  -hpisubt amp');  % for final MF, unless opted out of movecomp using mvcomp_fail above
cstr.tSSSstr = ' -st 10 -corr 0.98'; %'tSSSstr = '';
cstr.transtr = '';  % file for -trans, will be defined later

%% set up
% get subjects and runs to be processed (if specified)
if ~isempty(datinfo.subj_list)
    for ii = 1:length(datinfo.subj_list)
        sub_dirs_new{ii} = datinfo.sub_dirs{datinfo.subj_list(ii)};
        bad_chans_new = {};
        if ~isempty(datinfo.bad_chans)
            bad_chans_new{ii} = datinfo.bad_chans{datinfo.subj_list(ii)};
        end
        if isempty(datinfo.run_list)
            raw_files_in_new{ii} = datinfo.raw_files_in{datinfo.subj_list(ii)};
            raw_files_out_new{ii} = datinfo.raw_files_out{datinfo.subj_list(ii)};
        elseif isempty(datinfo.run_list{ii})
            raw_files_in_new{ii} = datinfo.raw_files_in{datinfo.subj_list(ii)};
            raw_files_out_new{ii} = datinfo.raw_files_out{datinfo.subj_list(ii)};
        else
            for jj = 1:length(datinfo.run_list{ii})
                raw_files_in_new{ii}{jj} = datinfo.raw_files_in{datinfo.subj_list(ii)}{datinfo.run_list{ii}(jj)};
                raw_files_out_new{ii}{jj} = datinfo.raw_files_out{datinfo.subj_list(ii)}{datinfo.run_list{ii}(jj)};
            end
        end
    end
    datinfo.sub_dirs = sub_dirs_new;
    datinfo.raw_files_in = raw_files_in_new;
    datinfo.raw_files_out = raw_files_out_new;
    datinfo.bad_chans = bad_chans_new;
end

datinfo.nr_subs = length(datinfo.sub_dirs);  % number of subdirectories/subjects


% Set up directory structures (only needs to be done once) and check if input files exist
file_flag = 0;
for ss = 1:datinfo.nr_subs  % create output directories if they don't exist yet
    now_dir = fullfile(datinfo.data_out_dir,datinfo.sub_dirs{ss});
    if exist(now_dir) ~= 7
        eval(sprintf('!mkdir %s -p', now_dir));
    end
    for rr = 1:length(datinfo.raw_files_in{ss})  % check if input fiff-files exist
        now_file  = fullfile(datinfo.data_in_dir, datinfo.sub_dirs{ss}, [datinfo.raw_files_in{ss}{rr} '.fif']);        
        if ~exist(now_file, 'file')
            fprintf(1, 'This input fiff-file does not exist: %s\n\n', now_file);
            file_flag = 1;
        end
    end
end
if file_flag  % don't continue if some files don't exist
    return
end

if ~isempty(datinfo.bad_chans) && (length(datinfo.bad_chans) ~= length(datinfo.sub_dirs))
    fprintf(1, '\nNumbers of bad channel lists (%d) and subject info (%d) do not correspond!\n', length(datinfo.bad_chans), length(datinfo.sub_dirs));
    return;
end;

if (length(datinfo.raw_files_in) ~= length(datinfo.raw_files_out)) && (length(datinfo.raw_files_out) ~= 1)
    fprintf(1, '\nNumbers of input (%d) and output (%d) lists do not match!\n', length(datinfo.raw_files_in), length(datinfo.raw_files_out));
    return;
end
%return

%% open matlabpool if required
% matlabpool close force CBU_Cluster
if ParType && datinfo.do_eval
    if matlabpool('size')==0
        MaxNsubs = 1;
        if ParType == 2
            for g=1:length(datinfo.nr_subs)
                MaxNsubs = max([MaxNsubs length(datinfo.sub_dirs{g})]);
            end
        end
        P = cbupool(MaxNsubs);
        P.SubmitArguments='-W x=\"NODESET:ONEOF:FEATURES:MAXFILTER\"';
        P.ResourceTemplate='-l nodes=^N^,mem=96GB,walltime=48:00:00';
        matlabpool(P);
    end
    fprintf(1, '\nMatlabpool initialised\n');
end

%% Main loop (can comment/uncomment lines below if want to parallelise over expts rather than subjects)
% display parameters just FYI
display( datinfo );
display( cstr);

if ParType == 0 || ParType == 1  % run non-parallel
    for s = 1:datinfo.nr_subs  % for all subjects

        maxfilter_one_subject(cstr, datinfo, s, ParType);

    end % subs
elseif ParType == 2,  % run parallel using parfor
    parfor s = 1:datinfo.nr_subs  % for all subjects

        maxfilter_one_subject(cstr, datinfo, s, ParType);

    end % subs
else,
    fprintf(1, 'Nothing to do for this ParType option\n');
end
    
%% Fin
if ParType
    matlabpool close force CBU_Cluster
end