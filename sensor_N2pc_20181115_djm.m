% N2pc
dbstop if error

% addpath(genpath('/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp'))
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
%addpath(genpath('/imaging/local/software/spm_toolbox/eeglab13_4_3b'))
spm('defaults', 'eeg');

workingdir = '/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp';

% Define SUBJECT INFORMATION
subs = [1,2,3,4,5,6,7,8,9,10,11,13,15,16,17,18,19,20];  % subject numbers
subjects_dirs = {'meg16_0317/161107','meg16_0319/161110','meg16_0321/161111','meg16_0322/161114','meg16_0325/161115','meg16_0327/161117','meg16_0330/161121','meg16_0332/161122','meg16_0333/161124','meg16_0337/161128','meg16_0339/161129','meg16_0340/161129','meg16_0341/161201','meg16_0343/161202','meg16_0345/161206','meg16_0346/161206','meg16_0348/161208','meg16_0349/161208','meg16_0350/161212','meg16_0352/161213'};
subjnum = [1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,6,2]; % counterbalancing numbers

modalities = {'EEG','MEG'};

group_N2pc_po7_con={[],[]};
group_N2pc_po7_ips={[],[]};
group_N2pc_po8_con={[],[]};
group_N2pc_po8_ips={[],[]};        
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
        
    for m = 1:2
        modality = modalities{m};
        if isequal(modality,'EEG')
            chan_po7 = indchannel(D,{'EEG066','EEG067','EEG071'  ,'EEG001','EEG052'});% PO7=66, PO3=67, O1=71; PO9=1, P7=52
            chan_po8 = indchannel(D,{'EEG070','EEG069','EEG073'  ,'EEG003','EEG060'});% PO8=70, PO4=69, O2=73; PO10=3, P8=60
        elseif isequal(modality,'MEG')
            chan_po7 = indchannel(D,{'MEG2142','MEG1933','MEG1922','MEG2043',...
                'MEG1742','MEG1733','MEG1942','MEG1913',...
                'MEG1712','MEG1723','MEG1642'});
            chan_po8 = indchannel(D,{'MEG2132','MEG2333','MEG2342','MEG2033',...
                'MEG2542','MEG2513','MEG2322','MEG2313',...
                'MEG2532','MEG2523','MEG2432'});
        end
        
        if sum(ismember(chan_po7,D.badchannels))<numel(chan_po7) ...
                && sum(ismember(chan_po8,D.badchannels))<numel(chan_po8)
            
            new_chan_po7 = setdiff(chan_po7,intersect(chan_po7,D.badchannels));
            
            N2pc_po7_con = (D.fttimelock.trial(TNiNc==-2,new_chan_po7,:));
            N2pc_po7_ips = (D.fttimelock.trial(TNiNc==-1,new_chan_po7,:));
            
            group_N2pc_po7_con{m}(end+1,:) = mean(mean(N2pc_po7_con,1),2);
            group_N2pc_po7_ips{m}(end+1,:) = mean(mean(N2pc_po7_ips,1),2);
            
            new_chan_po8 = setdiff(chan_po8,intersect(chan_po8,D.badchannels));
            
            N2pc_po8_con = (D.fttimelock.trial(TNiNc==-1,new_chan_po8,:));
            N2pc_po8_ips = (D.fttimelock.trial(TNiNc==-2,new_chan_po8,:));
            
            group_N2pc_po8_con{m}(end+1,:) = mean(mean(N2pc_po8_con,1),2);
            group_N2pc_po8_ips{m}(end+1,:) = mean(mean(N2pc_po8_ips,1),2);            
        end
        
        % %     po7 = D.selectdata('EEG066', [-0.1 1.5], {'image1'});
        % %     po8 = D.selectdata('EEG070', [-0.1 1.5], {'image1'});
    end % next modality
end %next subject
    
for m = 1:2    
    n=size(group_N2pc_po7_ips{m},1);
    % contralateral and ipsilateral
    figure(10+m); set(10+m,'position',[0,0,1420,950]); clf
    subplot(2,1,1)
    hold on;
    x1 = linspace(-100,1500,401);
    y1 = squeeze(mean(group_N2pc_po7_con{m})); % contra, left
    errd = squeeze(std(group_N2pc_po7_con{m}-group_N2pc_po7_ips{m},[],1))/sqrt(n);
    x1 = x1(1:150);
    y1 = y1(1:150);
    errd = errd(1:150);
    ph1=patch([x1 fliplr(x1)],[y1+errd fliplr(y1-errd)],[0 0 0.8],'facealpha',0.3,'edgecolor','none');
    plot(x1,y1,'LineWidth',1,'color',[0 0 0.8],'LineWidth',2);
    x2 = linspace(-100,1500,401);
    y2 = squeeze(mean(group_N2pc_po7_ips{m})); % ipsi, left
    x2 = x2(1:150);
    y2 = y2(1:150);
    ph2=patch([x2 fliplr(x2)],[y2+errd fliplr(y2-errd)],[0.8 0 0],'facealpha',0.3,'edgecolor','none');
    plot(x2,y2,'LineWidth',1,'color',[0.8 0 0],'LineWidth',2);
    axis tight
    plot([0,0],ylim,'k');
    plot(xlim,[0,0],'k');
    title('left hemisphere');

    subplot(2,1,2)
    hold on;
    x1 = linspace(-100,1500,401);
    y1 = squeeze(mean(group_N2pc_po8_con{m})); % contra, right
    errd = squeeze(std(group_N2pc_po8_con{m}-group_N2pc_po8_ips{m},[],1))/sqrt(n);
    x1 = x1(1:150);
    y1 = y1(1:150);
    errd = errd(1:150);
    ph1=patch([x1 fliplr(x1)],[y1+errd fliplr(y1-errd)],[0 0 0.8],'facealpha',0.3,'edgecolor','none');
    plot(x1,y1,'LineWidth',1,'color',[0 0 0.8],'LineWidth',2);
    x2 = linspace(-100,1500,401);
    y2 = squeeze(mean(group_N2pc_po8_ips{m})); % ipsi, right
    x2 = x2(1:150);
    y2 = y2(1:150);
    ph2=patch([x2 fliplr(x2)],[y2+errd fliplr(y2-errd)],[0.8 0 0],'facealpha',0.3,'edgecolor','none');
    plot(x2,y2,'LineWidth',1,'color',[0.8 0 0],'LineWidth',2);
    axis tight
    plot([0,0],ylim,'k');
    plot(xlim,[0,0],'k');
    title('right hemisphere');
    
    figure(20+m); clf
    hold on;
    x1 = linspace(-100,1500,401);
    y1 = (squeeze(mean(group_N2pc_po7_con{m})) + squeeze(mean(group_N2pc_po8_con{m})))/2;
    m1=(group_N2pc_po7_con{m} + group_N2pc_po8_con{m})/2;
    m2=(group_N2pc_po7_ips{m} + group_N2pc_po8_ips{m})/2;
    errd = squeeze(std(m1-m2,[],1))/sqrt(n);
    errd = errd(1:150);
    x1 = x1(1:150);
    y1 = y1(1:150);
    errd = errd(1:150);
    x2 = linspace(-100,1500,401);
    y2 = (squeeze(mean(group_N2pc_po7_ips{m})) + squeeze(mean(group_N2pc_po8_ips{m})))/2;
    x2 = x2(1:150);
    y2 = y2(1:150);
    xlim([-100,500])
    ylim([min(min(y1-errd,y2-errd)),max(max(y1+errd,y2+errd))])
    y_lim = ylim;
    patch('Faces',[1,2,3,4],'Vertices',[200,y_lim(1); 260,y_lim(1); 260,y_lim(2); 200,y_lim(2)],'FaceColor',[0.6,0.6,0.6],'FaceAlpha',0.4,'EdgeColor','none');
    
    addpath /imaging/dm01/MoreTools
    p1 = boundedline(x1,y1,errd,'cmap',[0 0 0.8],'alpha','transparency',0.2);
    p2 = boundedline(x2,y2,errd,'cmap',[0.8 0 0],'alpha','transparency',0.2);
    axis tight
%     p1 = plot(x1,y1,'LineWidth',1,'color',[0 0 0.8],'LineWidth',2);
%     p2 = plot(x2,y2,'LineWidth',1,'color',[0.8 0 0],'LineWidth',2);
%     ph1=patch([x1 fliplr(x1)],[y1+errd fliplr(y1-errd)],[0 0 0.8],'facealpha',0.2,'edgecolor','none');
%     ph2=patch([x2 fliplr(x2)],[y2+errd fliplr(y2-errd)],[0.8 0 0],'facealpha',0.2,'edgecolor','none');
%     title(sprintf('%s, both hemispheres combined',modalities{m}));
    legend([p1 p2], {'Contralateral','Ipsilateral'})

    plot([0,0],ylim,'k');
    plot(xlim,[0,0],'k');
    switch modalities{m}
        case 'EEG'
            ylabel('\muV')
        case 'MEG'
            ylabel('fT/cm^2')
    end
    xlabel('Time from stimulus onset (ms)')
    set(gca,'fontsize',20)
    
%     % FDR correction
%     [th, tp]=ttest(m1(:,1:150),m2(:,1:150));
%     addpath(genpath('/imaging/tw05/software'))
%     fdr=fdr_bh(tp);
%     plot(x1(logical(fdr)),zeros(size(find(fdr))),'.','color',[0 0 0],'MarkerSize',10);
    
%     % cluster correction
%     addpath(genpath('/imaging/tw05/Example Scripts/Myers et al/toolbox_plotting_and_stats'));
%     dat = m1(:,1:150)-m2(:,1:150);
%     nSims = 1000;
%     p_crit = 0.01;
%     p_thresh = 0.01;
%     [p,praw] = ClusterCorrection1(dat, nSims, p_crit, p_thresh);
%     sig_clusters = p<0.05;
%     plot(x1(logical(sig_clusters)),zeros(size(find(sig_clusters))),'.','color',[0 0 0],'MarkerSize',10);

      % TFCE
      addpath(genpath('/imaging/dm01/MoreTools/MatlabTFCE/'))
      [dp1, dp2]=matlab_tfce('onesample',2,permute(m1(:,1:150)-m2(:,1:150),[2,3,4,1]),[],[],1000,2,2/3,18,[],[]);
      sig=(dp1<0.025|dp2<0.025);
      plot(x1(logical(sig)),zeros(size(find(sig))),'.','color',[0 0 0],'MarkerSize',10);


      cd('/imaging/tw05/Preparatory_Attention_Study/Version3-FullExp/group_results')
      saveas(gcf,sprintf('posterior_sensor_N2pc_%s.png',modalities{m}))
      print(figure(20+m),sprintf('posterior_sensor_N2pc_%s.eps',modalities{m}),'-depsc2','-painters');
end % EEG then MEG

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



%% t-test

% po7_hvector = zeros(1,401); % hypothesis
% po7_pvector = zeros(1,401); % p-values
% po8_hvector = zeros(1,401); % hypothesis
% po8_pvector = zeros(1,401); % p-values
% for time = 1:401
%     [po7_hvector(time) po7_pvector(time)] = ttest2(squeeze(group_N2pc_po7_con(:,time)),squeeze(group_N2pc_po7_ips(:,time)));
%     [po8_hvector(time) po8_pvector(time)] = ttest2(squeeze(group_N2pc_po8_con(:,time)),squeeze(group_N2pc_po8_ips(:,time)));
% end
%
%
% % fdr correction
% addpath(genpath('/imaging/tw05/Preparatory_Attention_Study/software'))
% fdr(1,:)=fdr_bh(pvals(1,:),0.05);
% fdr(2,:)=fdr_bh(pvals(2,:),0.05);
% fdr(3,:)=fdr_bh(pvals(3,:),0.05);
% fdr=double(fdr);
% fdr(fdr==0) = -999;
