
function maxfilter_one_subject(cstr, datinfo, s, ParType)
% run Maxfilter on sessions within sub-directory for one subject
% cstr: structure with maxfilter commands
% datinfo: data information (files, directories etc.)
% s: subject index to process here
% ParType: Parallelisation option
% Olaf Hauk, Sep 14 (adapted from RH script)

% for ease of use, get info from datinfo structure                      
data_in_dir = datinfo.data_in_dir;
data_out_dir = datinfo.data_out_dir;
sub_dirs = datinfo.sub_dirs;
nr_subs = datinfo.nr_subs;
bad_chans = datinfo.bad_chans;                
OverWrite = datinfo.OverWrite;
movfile = datinfo.movfile;        
Nbuf = datinfo.Nbuf;           

if length(datinfo.raw_files_in) > 1  % if different input file names provided for each subject
    raw_files_in = datinfo.raw_files_in{s};
else  % if only one list of file names provided for all subjects
    raw_files_in = datinfo.raw_files_in{1};
end
if length(datinfo.raw_files_out) > 1  % if different output file names provided for each subject
    raw_files_out = datinfo.raw_files_out{s};
else  % if only one list of file names provided for all subjects
    raw_files_out = datinfo.raw_files_out{1};
end
nr_raw_files = length(raw_files_in);

if ~isempty(datinfo.mvcomp_fail)  % if individual movecomp failures specified
    mvcomp_fail = datinfo.mvcomp_fail{s};
else  % if no movecomp failures specified
    mvcomp_fail = zeros(1,length(raw_files_in));
end;

finstr = [];  % initialise final Maxfilter command string

% Output directory
out_dir = fullfile(data_out_dir, sub_dirs{s}), 

if datinfo.do_eval  % delete/create file for movement parameters
    if exist(movfile)
        delete(movfile)
    end
    eval(sprintf('!touch %s',movfile));            
end

log_filename = [num2str(s) '_' datinfo.logfilestem];
log_fid = fopen( log_filename, 'w' );

fprintf(log_fid, '%s\n', sub_dirs{s});

for r = 1:nr_raw_files            

   raw_file = fullfile(data_in_dir,sub_dirs{s},sprintf('/%s.fif',raw_files_in{r}));  % Get raw FIF file

   fprintf(log_fid, '\n%s\n', raw_file);

   raw_stem = raw_file(1:(end-4));

   badstr = cstr.badstr;    % will be altered during loop

   % determine sphere origin
   % Fit sphere (since better than MaxFilter does)
   if isempty(cstr.headori)
      if r == 1  % fit sphere doesn't change with run!
            incEEG = 0;
%                     if exist(fullfile(out_dir,'fittmp.txt')); delete(fullfile(out_dir,'fittmp.txt')); end
%                     if exist(fullfile(out_dir,sprintf('run_%02d_hpi.txt',raw_files_out{r}))); delete(fullfile(out_dir,sprintf('run_%02d_hpi.txt',raw_files_out{r})));  end
%                     [orig,rad,fit] = meg_fit_sphere(raw_file,out_dir,sprintf('%s_hpi.txt',raw_files_out{r}),incEEG);
%                     delete(fullfile(out_dir,'fittmp.txt'));
            if exist(fullfile(out_dir,sprintf('run_%02d_hpi.txt',raw_files_out{r}))); delete(fullfile(out_dir,sprintf('run_%02d_hpi.txt',raw_files_out{r})));  end
            cwd = pwd;
            [orig,rad,fit] = meg_fit_sphere(raw_file,out_dir,sprintf('%s_hpi.txt',raw_files_out{r}),incEEG);
            cd( cwd );  % back to original working directory
        end
        cstr.origstr = sprintf(' -origin %d %d %d -frame head',orig(1),orig(2),orig(3));
    else,
        % use pre-defined origin
        cstr.origstr = sprintf(' -origin %d %d %d -frame head',cstr.headori(1),cstr.headori(2),cstr.headori(3))
    end;

    fprintf(1, '\nHead origin: %s\n', cstr.origstr);


    %% 1. Bad channel detection (this email says important if doing tSSS later https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=NEUROMEG;d3f363f3.1205)

    outfile = fullfile(out_dir,sprintf('%s_bad',raw_files_out{r}));
    badfile = sprintf('%s.txt',outfile);
    bad_filename = [outfile '.fif'];
    Dout = dir( bad_filename );
    % if output file exists, get date
    if isempty( Dout )
        bad_datenum = 0;
    else
        bad_datenum = datenum( Dout.date )
    end

    cstr.filestr = sprintf(' -f %s -o %s.fif',raw_file,outfile);

    finstr_bad = [cstr.maxfstr cstr.filestr cstr.origstr cstr.basestr cstr.badstr cstr.compstr_bad];

    % get bad channels from maxfilter
    bad_str_tmp = maxfilter_get_bad_channels_from_fiff(raw_file, raw_files_out{r}, badfile, out_dir, finstr_bad, Nbuf, OverWrite, ParType, datinfo.do_eval); 



    % Take note whether output files were successfully created or not
    Dout = dir( bad_filename );
    if isempty( Dout )
        bad_datenum2 = 0;
    else
        bad_datenum2 = datenum( Dout.date )
    end
    if (bad_datenum2 > bad_datenum),
        fprintf(log_fid, 'Created: %s\n', bad_filename);
        fprintf(log_fid, '%d  %s\n', Dout.bytes, Dout.date);
    else
        fprintf(log_fid, 'NOT created: %s\n', bad_filename);
    end

    % delete current fiff-file if desired
    if datinfo.del_bads && exist(bad_filename, 'file')  
        delete( bad_filename );
    end

    % add manually specified bad channels
    if ~isempty(bad_chans)
        if isempty(bad_str_tmp) && ~isempty(bad_chans{s}),  % if -bad wasn't added during bad channel detection
            bad_str_tmp = ' -bad ';
        end            
        badstr = [badstr bad_str_tmp];
        for bb = 1:length(bad_chans{s}),
            badstr = [badstr ' ' bad_chans{s}{bb} ' '];                            
        end
    end
    fprintf(1, 'Bad channels: %s\n\n', badstr);                        

    %% 2. tSSS (incl. align within subject if multiple runs)

    outfile = fullfile(out_dir,sprintf('%s', raw_files_out{r}));

    % trans to first file in list
    if r == 1
        transfstfile = [outfile '.fif'];
        cstr.transtr = '';
    else
        cstr.transtr = sprintf(' -trans %s ', transfstfile);
    end

    % get date of existing file
    Dout = dir( [outfile '.fif'] );
    if isempty( Dout )
        out_datenum = 0;
    else
        out_datenum = datenum( Dout.date )
    end            

    if ~exist(sprintf('%s.fif',outfile)) | OverWrite                            

        cstr.filestr = sprintf(' -f %s -o %s.fif', raw_file,outfile);                                                                               

        % check if movement compensation shall be applied
        if mvcomp_fail, % if not
            finstr = [cstr.maxfstr cstr.filestr cstr.basestr badstr cstr.tSSSstr cstr.origstr cstr.transtr cstr.dsstr sprintf(' -v | tee %s.log',outfile)]
            fprintf(1, 'No movement compensation applied for %s\n', outfile);
            fprintf(log_fid, 'No movement compensation requested for %s\n', outfile);
        else,  % if yes
            finstr = [cstr.maxfstr cstr.filestr cstr.basestr badstr cstr.tSSSstr cstr.compstr cstr.origstr cstr.transtr sprintf(' -v | tee %s.log',outfile)]
        end; 

        fprintf(1, '\nFinal Maxfilter command:\n%s\n\n', finstr);
        if r == 1
            fprintf(log_fid, '%s\n\n', finstr);
        end;
        
        if datinfo.do_eval  % only execute if desired
            if ParType==1
                spmd; eval(finstr); end
            else
                eval(finstr);
            end                

            eval(sprintf('!echo ''Trans 1st...'' >> %s', movfile));
            eval(sprintf('!cat %s.log | sed -n ''/Position change/p'' | cut -f 7- -d '' '' >> %s', outfile,movfile));
        end
    end

    % Take note whether output files were successfully created or not
    % try again without movecomp first
    Dout = dir( [outfile '.fif'] );
    if isempty( Dout )
        out_datenum2 = 0;
    else
        out_datenum2 = datenum( Dout.date )
    end  

    if datinfo.do_eval  % only execute if desired
        if (out_datenum2 == out_datenum) && OverWrite, % if date of output file has not changed when it should have...
            finstr = [cstr.maxfstr cstr.filestr cstr.basestr badstr cstr.tSSSstr cstr.origstr cstr.transtr cstr.dsstr sprintf(' -v | tee %s.log',outfile)]
            fprintf(1, 'Trying again without movement compensation for %s\n', outfile);
            fprintf(log_fid, 'NOT created at first attempt: %s\n', [outfile '.fif']);                
            if ParType==1
                spmd; eval(finstr); end;                                               
            else
                eval(finstr);
            end 

            Dout = dir( [outfile '.fif'] );
            if isempty( Dout )
                out_datenum2 = 0;
            else
                out_datenum2 = datenum( Dout.date )
            end 
            if out_datenum2 > out_datenum,
                fprintf(log_fid, 'Created (2nd attempt): %s\n', [outfile '.fif']);
                fprintf(log_fid, '%d  %s\n', Dout.bytes, Dout.date);                
            else
                fprintf(log_fid, 'NOT created (2nd attempt): %s\n', [outfile '.fif']);
            end
        elseif (out_datenum2 > out_datenum) && (out_datenum>0)
            fprintf(log_fid, 'Created and overwritten: %s\n', [outfile '.fif']);
        elseif (out_datenum2 > out_datenum) && (out_datenum==0)
            fprintf(log_fid, 'Created new: %s\n', [outfile '.fif']);
        elseif (out_datenum > 0) && ~OverWrite,
            fprintf(log_fid, 'Existing not overwritten: %s\n', [outfile '.fif']);
        elseif (out_datenum == 0) && (out_datenum2==0),
            frintf(log_fid, 'Nothing created for: %s\n', [outfile '.fif']);
        else
            fprintf(log_fid, 'Not sure what do to with: %s\n', [outfile '.fif']);
        end

    end % if datinfo...   
    
    % if downsampling specified
    if ~isempty(cstr.dsstr)
        outfile_ds = [outfile '_ds'];
        cstr.filestr = sprintf(' -f %s.fif -o %s.fif', outfile, outfile_ds);
        
        Dout = dir( [outfile_ds '.fif'] );
        if isempty( Dout )
            out_datenum = 0;
        else
            out_datenum = datenum( Dout.date )
        end 
        
        finstr = [cstr.maxfstr cstr.filestr cstr.dsstr sprintf(' -nosss -force  -v | tee %s.log',outfile_ds)];
        if datinfo.do_eval
            fprintf(1, '\nDownsampling to: %s\n', outfile_ds);
            eval(finstr);
        else
            fprintf(1, '\nDownsampling command:\n %s\n', finstr);
        end
        
        Dout = dir( [outfile_ds '.fif'] );
        if isempty( Dout )
            out_datenum2 = 0;
        else
            out_datenum2 = datenum( Dout.date )
        end 
        if out_datenum2 > out_datenum,
            fprintf(log_fid, 'Successfully downsampled: %s\n', [outfile_ds '.fif']);
            fprintf(log_fid, '%d  %s\n', Dout.bytes, Dout.date);                
        else
            fprintf(log_fid, 'NOT downsampled: %s\n', [outfile_ds '.fif']);
        end
    end

end  % for r...

fprintf(log_fid, '\n\nI''ve done everything I could!');

fclose( log_fid );