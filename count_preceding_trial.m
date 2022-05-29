% attempt to solve TNiNc baseline decoding
clear
load('CrossGenMatrix032msRaw_sub1.mat');
goodtrials = indtrial(D,D.condlist,'GOOD');

% things that proceed T
T_ind = intersect(goodtrials,find(TNiNc==1)) ;
before_T_ind = TNiNc(T_ind-1);
before_T_T = sum(before_T_ind==1);
before_T_Ni = sum(before_T_ind==2);
before_T_Nc = sum(before_T_ind==3);
before_T_3itemCue1 = sum(before_T_ind==-1);
before_T_3itemCue2 = sum(before_T_ind==-2);
before_T_bad = sum(before_T_ind==0);
% things that proceed Ni
Ni_ind = intersect(goodtrials,find(TNiNc==2));
before_Ni_ind = TNiNc(Ni_ind-1);
before_Ni_T = sum(before_Ni_ind==1);
before_Ni_Ni = sum(before_Ni_ind==2);
before_Ni_Nc = sum(before_Ni_ind==3);
before_Ni_3itemCue1 = sum(before_Ni_ind==-1);
before_Ni_3itemCue2 = sum(before_Ni_ind==-2);
before_Ni_bad = sum(before_Ni_ind==0);

