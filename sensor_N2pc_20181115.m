% N2pc
dbstop if error

% addpath(genpath('/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp'))
% addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'))
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/'))
addpath(genpath('/imaging/local/software/spm_toolbox/eeglab13_4_3b'))
spm('defaults', 'eeg');

workingdir = '/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp';

% Define SUBJECT INFORMATION
subs = [1,2,3,4,5,6,7,8,9,10,11,13,15,16,17,18,19,20];  % subject numbers
subjects_dirs = {'meg16_0317/161107','meg16_0319/161110','meg16_0321/161111','meg16_0322/161114','meg16_0325/161115','meg16_0327/161117','meg16_0330/161121','meg16_0332/161122','meg16_0333/161124','meg16_0337/161128','meg16_0339/161129','meg16_0340/161129','meg16_0341/161201','meg16_0343/161202','meg16_0345/161206','meg16_0346/161206','meg16_0348/161208','meg16_0349/161208','meg16_0350/161212','meg16_0352/161213'};
subjnum = [1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,6,2]; % counterbalancing numbers

modalities = {'EEG','MEG'};
for m = 1:2
    modality = modalities{m};
    subind_po7 = 1; subind_po8 = 1;
    for sub = subs
        
        cd(workingdir)
        swd = sprintf('sub%02d/%s',sub,subjects_dirs{sub}); % subject working directory
        cd(swd)
        D = spm_eeg_load('caefMattn2_attention_task_block1_raw.mat');
        
        
        % load the counterbalacing
        % row indicates the 6 possibilities
        % column indicates the position of (col:1)face, (col:2)house, (col:3)violin;
        % 1 is top-left, 2 is top-right, 3 is bottom-middle
        pos_all = perms(1:3);
        
        % define subject specific pairings
        allperms = perms([1 2 3]);
        for i = 1:length(allperms)
            if rem(subjnum(sub),6)+1 == i
                choosetargets = allperms(i,:);
            end
        end
        target1 = sprintf('image%d',choosetargets(1));
        target2 = sprintf('image%d',choosetargets(2));
        target3 = sprintf('image%d',choosetargets(3));
        
        t1_left = find(pos_all(:,choosetargets(1))==1); % top-left == 1
        t1_right = find(pos_all(:,choosetargets(1))==2); % top-right == 2
        t2_left = find(pos_all(:,choosetargets(2))==1); % top-left == 1
        t2_right = find(pos_all(:,choosetargets(2))==2); % top-right == 2
        
        % load condition list
        conds = D.conditions;
        
        % calculate list of 1-items and 3-items
        TNiNc = zeros(1,length(conds));
        indexcue1 = strfind(conds,'cue1');
        indexcue1 = ~cellfun(@isempty,indexcue1);
        for findcue1 = 1:length(indexcue1)
            if indexcue1(findcue1) == 1
                for i = 1:3 %three events per trial
                    try % use try because if the last cue does not follow three trials it will crash
                        if strcmp(conds(findcue1+i),'cue1')==1 || strcmp(conds(findcue1+i),'cue2')==1
                            break
                        elseif strcmp(conds(findcue1+i),target1) == 1
                            TNiNc(findcue1+i) = 1; %target
                        elseif strcmp(conds(findcue1+i),target2) == 1
                            TNiNc(findcue1+i) = 2; %inconsistent non-target
                        elseif strcmp(conds(findcue1+i),target3) == 1
                            TNiNc(findcue1+i) = 3; %consistent non-target
                        elseif regexp(conds{findcue1+i},sprintf('imageall_pos%d',t1_left(1))) == 1 | regexp(conds{findcue1+i},sprintf('imageall_pos%d',t1_left(2))) == 1
                            TNiNc(findcue1+i) = -1; %three-item display (target left visual field)
                        elseif regexp(conds{findcue1+i},sprintf('imageall_pos%d',t1_right(1))) == 1 | regexp(conds{findcue1+i},sprintf('imageall_pos%d',t1_right(2))) == 1
                            TNiNc(findcue1+i) = -2; %three-item display (target right visual field)
                        else TNiNc(findcue1+i) = 0; %opaque
                        end
                    catch
                    end
                end
            end
        end
        indexcue2 = strfind(conds,'cue2');
        indexcue2 = ~cellfun(@isempty,indexcue2);
        for findcue2 = 1:length(indexcue2)
            if indexcue2(findcue2) == 1
                for i = 1:3 %three events per trial
                    try
                        if strcmp(conds(findcue2+i),'cue1')==1 || strcmp(conds(findcue2+i),'cue2')==1
                            break
                        elseif strcmp(conds(findcue2+i),target2) == 1
                            TNiNc(findcue2+i) = 1; %target
                        elseif strcmp(conds(findcue2+i),target1) == 1
                            TNiNc(findcue2+i) = 2; %inconsistent non-target
                        elseif strcmp(conds(findcue2+i),target3) == 1
                            TNiNc(findcue2+i) = 3; %consistent non-target
                        elseif regexp(conds{findcue2+i},sprintf('imageall_pos%d',t2_left(1))) == 1 | regexp(conds{findcue2+i},sprintf('imageall_pos%d',t2_left(2))) == 1
                            TNiNc(findcue2+i) = -1; %three-item display (target left visual field)
                        elseif regexp(conds{findcue2+i},sprintf('imageall_pos%d',t2_right(1))) == 1 | regexp(conds{findcue2+i},sprintf('imageall_pos%d',t2_right(2))) == 1
                            TNiNc(findcue2+i) = -2; %three-item display (target right visual field)
                        else TNiNc(findcue2+i) = 0; %opaque
                        end
                    catch
                    end
                end
            end
        end
        
        bad_trials = D.badtrials;
        TNiNc(bad_trials) = 99;
        
        if isequal(modality,'EEG')
            chan_po7 = indchannel(D,{'EEG066','EEG067','EEG071'});
            chan_po8 = indchannel(D,{'EEG070','EEG069','EEG073'});
        elseif isequal(modality,'MEG')
            chan_po7 = indchannel(D,{'MEG2142','MEG1933','MEG1922','MEG2043',...
                'MEG1742','MEG1733','MEG1942','MEG1913',...
                'MEG1712','MEG1723','MEG1642'});
            chan_po8 = indchannel(D,{'MEG2132','MEG2333','MEG2342','MEG2033',...
                'MEG2542','MEG2513','MEG2322','MEG2313',...
                'MEG2532','MEG2523','MEG2432'});
        end
        
        if sum(ismember(chan_po7,D.badchannels))<numel(chan_po7)
            
            new_chan_po7 = setdiff(chan_po7,intersect(chan_po7,D.badchannels));
            
            N2pc_po7_con = D.fttimelock.trial(TNiNc==-2,new_chan_po7,:);
            N2pc_po7_ips = D.fttimelock.trial(TNiNc==-1,new_chan_po7,:);
            
            group_N2pc_po7_con(subind_po7,:) = mean(mean(N2pc_po7_con,1),2);
            group_N2pc_po7_ips(subind_po7,:) = mean(mean(N2pc_po7_ips,1),2);
            
            subind_po7 = subind_po7 + 1;
            
        end
        
        if sum(ismember(chan_po8,D.badchannels))<numel(chan_po8)
            
            new_chan_po8 = setdiff(chan_po8,intersect(chan_po8,D.badchannels));
            
            N2pc_po8_con = D.fttimelock.trial(TNiNc==-1,new_chan_po8,:);
            N2pc_po8_ips = D.fttimelock.trial(TNiNc==-2,new_chan_po8,:);
            
            group_N2pc_po8_con(subind_po8,:) = mean(mean(N2pc_po8_con,1),2);
            group_N2pc_po8_ips(subind_po8,:) = mean(mean(N2pc_po8_ips,1),2);
            
            subind_po8 = subind_po8 + 1;
            
        end
        
        % %     po7 = D.selectdata('EEG066', [-0.1 1.5], {'image1'});
        % %     po8 = D.selectdata('EEG070', [-0.1 1.5], {'image1'});
    end
    
    cd /imaging/tw05/Preparatory_Attention_Study/Version3-FullExp/paper_figures
    % contralateral and ipsilateral
    figure
    set(gca,'fontsize',12);
    subplot(2,1,1)
    hold on;
    x1 = linspace(-100,1500,401);
    y1 = squeeze(mean(group_N2pc_po7_con));
    err1 = squeeze(std(group_N2pc_po7_con,[],1))/sqrt(subind_po7);
    x1 = x1(1:150);
    y1 = y1(1:150);
    err1 = err1(1:150);
    %patch([x1 fliplr(x1)],[y1+err1 fliplr(y1-err1)],[0.3 0 0],'facealpha',0.4,'edgecolor','none');
    plot(x1,y1,'LineWidth',1,'color',[1 0 0],'LineWidth',2);
    x2 = linspace(-100,1500,401);
    y2 = squeeze(mean(group_N2pc_po7_ips));
    err2 = squeeze(std(group_N2pc_po7_ips,[],1))/sqrt(subind_po7);
    x2 = x2(1:150);
    y2 = y2(1:150);
    err2 = err2(1:150);
    %patch([x2 fliplr(x2)],[y2+err2 fliplr(y2-err2)],[0 0 0.3],'facealpha',0.4,'edgecolor','none');
    plot(x2,y2,'LineWidth',1,'color',[0 0 1],'LineWidth',2);
    plot([0,0],[min(y1)-2,max(y1)+2],'k');
    plot([min(x1),max(x1)],[0,0],'k');
    ylim([min(y1)-2,max(y1)+2]);
    title('left hemisphere');
    subplot(2,1,2)
    hold on;
    x1 = linspace(-100,1500,401);
    y1 = squeeze(mean(group_N2pc_po8_con));
    err1 = squeeze(std(group_N2pc_po8_con,[],1))/sqrt(subind_po8);
    x1 = x1(1:150);
    y1 = y1(1:150);
    err1 = err1(1:150);
    %patch([x1 fliplr(x1)],[y1+err1 fliplr(y1-err1)],[0.3 0 0],'facealpha',0.4,'edgecolor','none');
    plot(x1,y1,'LineWidth',1,'color',[1 0 0],'LineWidth',2);
    x2 = linspace(-100,1500,401);
    y2 = squeeze(mean(group_N2pc_po8_ips));
    err2 = squeeze(std(group_N2pc_po8_ips,[],1))/sqrt(subind_po8);
    x2 = x2(1:150);
    y2 = y2(1:150);
    err2 = err2(1:150);
    %patch([x2 fliplr(x2)],[y2+err2 fliplr(y2-err2)],[0 0 0.3],'facealpha',0.4,'edgecolor','none');
    plot(x2,y2,'LineWidth',1,'color',[0 0 1],'LineWidth',2);
    plot([0,0],[min(y1)-2,max(y1)+2],'k');
    plot([min(x1),max(x1)],[0,0],'k');
    ylim([min(y1)-2,max(y1)+2]);
    title('right hemisphere');
    saveas(gcf,sprintf('sensor_N2pc_unilateral_%s.png',modality))
    
    figure
    set(gca,'fontsize',12);
    hold on;
    x1 = linspace(-100,1500,401);
    y1 = (squeeze(mean(group_N2pc_po7_con)) + squeeze(mean(group_N2pc_po8_con)))/2;
    err1 = squeeze(std((group_N2pc_po7_con + group_N2pc_po8_con)/2,[],1))/sqrt(subind_po7);
    x1 = x1(1:150);
    y1 = y1(1:150);
    err1 = err1(1:150);
    plot(x1,y1,'LineWidth',1,'color',[1 0 0],'LineWidth',2);
    x2 = linspace(-100,1500,401);
    y2 = (squeeze(mean(group_N2pc_po7_ips)) + squeeze(mean(group_N2pc_po8_ips)))/2;
    err2 = squeeze(std((group_N2pc_po7_ips + group_N2pc_po8_ips)/2,[],1))/sqrt(subind_po7);
    x2 = x2(1:150);
    y2 = y2(1:150);
    err2 = err2(1:150);
    plot(x2,y2,'LineWidth',1,'color',[0 0 1],'LineWidth',2);
    plot([0,0],[min(y1)-2,max(y1)+2],'k');
    plot([min(x1),max(x1)],[0,0],'k');
    if isequal(modality,'EEG')
    ylim([min(y1)-2,max(y1)+2]);
    else
    ylim([min(y1)-1,max(y1)+1]);
    end
    title('both hemispheres combined');
    saveas(gcf,sprintf('sensor_N2pc_bilateral_%s.png',modality))
end
% difference waves
% figure
% subplot(2,1,1)
% hold on;
% x1 = linspace(-100,1500,401);
% y1 = squeeze(mean(group_N2pc_po7_con))'-squeeze(mean(group_N2pc_po7_ips))';
% err1 = squeeze(std((group_N2pc_po7_con-group_N2pc_po7_ips),[],1))'/sqrt(numel(subs));
% patch([x1 fliplr(x1)],[y1+err1 fliplr(y1-err1)],[0 0.3 0],'facealpha',0.4,'edgecolor','none');
% plot(x1,y1,'LineWidth',1,'color',[0 1 0]);
% title('PO7');
% subplot(2,1,2)
% hold on;
% x1 = linspace(-100,1500,401);
% y1 = squeeze(mean(group_N2pc_po8_con))'-squeeze(mean(group_N2pc_po8_ips))';
% err1 = squeeze(std((group_N2pc_po8_con-group_N2pc_po8_ips),[],1))'/sqrt(numel(subs));
% patch([x1 fliplr(x1)],[y1+err1 fliplr(y1-err1)],[0 0.3 0],'facealpha',0.4,'edgecolor','none');
% plot(x1,y1,'LineWidth',1,'color',[0 1 0]);
% title('PO8');

cd('/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp/group_results')
save('posterior_sensor_N2pc.mat')
saveas(gcf,'posterior_sensor_N2pc.png')


%% t-test

po7_hvector = zeros(1,401); % hypothesis
po7_pvector = zeros(1,401); % p-values
po8_hvector = zeros(1,401); % hypothesis
po8_pvector = zeros(1,401); % p-values
for time = 1:401
    [po7_hvector(time) po7_pvector(time)] = ttest2(squeeze(group_N2pc_po7_con(:,time)),squeeze(group_N2pc_po7_ips(:,time)));
    [po8_hvector(time) po8_pvector(time)] = ttest2(squeeze(group_N2pc_po8_con(:,time)),squeeze(group_N2pc_po8_ips(:,time)));
end
%
%
% % fdr correction
% addpath(genpath('/imaging/tw05/Preparatory_Attention_Study/software'))
% fdr(1,:)=fdr_bh(pvals(1,:),0.05);
% fdr(2,:)=fdr_bh(pvals(2,:),0.05);
% fdr(3,:)=fdr_bh(pvals(3,:),0.05);
% fdr=double(fdr);
% fdr(fdr==0) = -999;
