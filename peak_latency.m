%% Cross-validation of which target the participant is attending for during the preparatory phase. Done separately within each ROI
% written by Tanya Wen
% 2016/10/13

dbstop if error

% addpath(genpath('/imaging/tw05'))
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'))
addpath(genpath('/imaging/local/software/spm_toolbox/eeglab13_4_3b'))
spm('defaults', 'eeg');

workingdir = '/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp';
% Define SUBJECT INFORMATION
% Define SUBJECT INFORMATION
subs = [1,2,3,4,5,6,7,8,9,10,11,13,15,16,17,18,19,20];  % subject numbers
subjects_dirs = {'meg16_0317/161107','meg16_0319/161110','meg16_0321/161111','meg16_0322/161114','meg16_0325/161115','meg16_0327/161117','meg16_0330/161121','meg16_0332/161122','meg16_0333/161124','meg16_0337/161128','meg16_0339/161129','meg16_0340/161129','meg16_0341/161201','meg16_0343/161202','meg16_0345/161206','meg16_0346/161206','meg16_0348/161208','meg16_0349/161208','meg16_0350/161212','meg16_0352/161213'};
subjnum = [1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,6,2]; % counterbalancing numbers

roi_list = {'MD_ESV.nii','MD_IPS.nii','MD_AI.nii','MD_postMFG.nii','MD_ACCpreSMA.nii',...
    'MD_FEF.nii','MD_antMFG.nii','MD_midMFG.nii','Auditory_Te3.nii','Visual_hOc1.nii'};
roinames=roi_list;
roinames = strrep(roinames, '_', '-');
roinames = strrep(roinames, '.nii', '');

subind=1;
for sub = subs
    
    % go to subject directory
    cd(workingdir)
    swd = sprintf('sub%02d/%s',sub,subjects_dirs{sub}); % subject working directory
    cd(fullfile(swd,'crossval_ROIs_all'))

    for rois = 1:numel(roi_list)
        load(sprintf('CrossVal032msRaw_%s_sub%s.mat',roinames{rois},num2str(sub)));
        accuracy_sub(subind,rois,:) = mean(accuracy_vals(:,:),1); 
    end
    
    subind = subind + 1;
end


cd(workingdir)
output_dir = 'group_results';
try cd(output_dir)
catch eval(sprintf('!mkdir %s',output_dir)); cd(output_dir);
end

window = [1 size(data3D,3)]; % time of interest
wins      = window(1):stpsize:window(2);

subs = [1,2,3,4,5,6,7,8,9,10,11,13,15,16,17,18,19,20];  % subject numbers
figure('position', [0, 0, 600, 600])  % create new figure with specified size
for sub = 1:length(subs)
    for rois = 1:numel(roi_list)
        
        trial_length = size(accuracy_vals,2);
        [value peak(sub,rois)] = max(accuracy_sub(sub,rois,1:trial_length/2));
        timing(sub,rois) = wins(peak(sub,rois))*4-100;
        
    end
end
bar(mean(timing,1));
hold on
errorbar(mean(timing,1),std(timing,[],1)/sqrt(numel(subs)),'.');
set(gca,'XTickLabel',roinames)
set(gca,'XTickLabelRotation',45)
ylabel('milliseconds')
title('peak decoding')
saveas(gcf,sprintf('Group-PeakLatency032ms.png'));
close all


% group peak
for rois = 1:numel(roi_list)
    y = squeeze(mean(accuracy_sub(:,rois,:),1))';
    [value peak(rois)] = max(y);
    timing(rois) = wins(peak(rois))*4-100;
end
