function badstr = bad_channels_from_mf(raw_file, raw_file_out, badfile, out_dir, finstr, Nbuf, OverWrite, ParType, do_eval)
% get bad channels for one fiff-file from maxfilter using -autobad
% raw_file: fiff-file for which bad channels needed
% raw_file_out: name of output fiff-file from maxfilter
% badfile: file that records bad channel information
% out_dir: output directory
% finstr: text string with maxfilter parameters
% Nbuf: number of buffers from fiff-file
% OverWrite: 0/1, whether output files to be overwritten
% ParType: Parallelisation option
% do_eval: 0/1, whether commands to be executed, or only to be displayed at command line
% Olaf Hauk, Sep 14 (adapated from RH script)

outfile = fullfile(out_dir,sprintf('%s_bad',raw_file_out));

if ~exist(badfile) | OverWrite                         
                       
    % Write out movements too...
    posfile = fullfile(out_dir,sprintf('%s_headpos.txt',raw_file_out));

    finstr = [finstr ' -hp ' posfile];

    finstr = [finstr ' -v | tee ' outfile '.log'];

    fprintf(1, '\n\nMF command for bad channel detection: \n%s\n\n', finstr);

    if do_eval
        if ParType==1
            spmd; eval(finstr); 
        end
        else
            eval(finstr);
        end
    end

    % Pull out bad channels from logfile:
    cat_bad_cmd = sprintf('!cat %s.log | sed -n -e ''/Detected/p'' -e ''/Static/p'' | cut -f 5- -d '' '' >> %s',outfile,badfile);
    fprintf(1, '\n\nNow attempting:\n%s\n\n', cat_bad_cmd);
    fp = fopen(badfile,'w'); fprintf(fp,' '); fclose(fp);
    if do_eval
        eval(cat_bad_cmd);
    end
end


fprintf(1, '\nDetermining bad channels\n');
tmp=dlmread(badfile,' '); % Nbuf = size(tmp,1);
tmp=reshape(tmp,1,prod(size(tmp)));
tmp=tmp(tmp>0); % Omit zeros (padded by dlmread):

% Get frequencies (number of buffers in which chan was bad):
[frq,allbad] = hist(tmp,unique(tmp));

% Mark bad based on threshold (currently ~5% of buffers (assuming 500 buffers)):
badchans = allbad(frq>0.05*Nbuf);
if isempty(badchans)
    badstr = '';
else
    badstr = sprintf(' -bad %s',num2str(badchans))
end

end
