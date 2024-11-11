
%% Figures and Statistics Script

load('BreathworkCO2-RawData.mat')
% you can also read in data directly from Excel using the
% BreathworkCO2_AllRawDataImports script

noftimepoints = 6;
sessionpeakwindow = 2:4;
nofpassive = length(find(breathvector==1)); %18
nofactive = length(find(breathvector==2)); %43
ntotal = length(breathvector);%61
nholo = length(find(holconvector==1)); %30
nconcon = length(find(holconvector==2)); %31

%% Figure 1a-b
co2means = nanmean(co2times(:,sessionpeakwindow)')';
co2minimums = min(co2times(:,sessionpeakwindow)')';

figure;
subplot(1,2,1)
boxplot(co2means, [breathvector,holconvector])
ylim([0 55])
subplot(1,2,2)
boxplot(co2minimums, [breathvector,holconvector])
ylim([0 55])

% nanmean(co2means(breathvector==1))
% nanstd(co2means(breathvector==1))/sqrt(nofpassive)
% nanmean(co2means(breathvector==2))
% nanstd(co2means(breathvector==2))/sqrt(nofactive)
% 
% nanmean(co2minimums(breathvector==1))
% nanstd(co2minimums(breathvector==1))/sqrt(nofpassive)
% nanmean(co2minimums(breathvector==2))
% nanstd(co2minimums(breathvector==2))/sqrt(nofactive)

[p,tbl,stats] = anovan(co2means,{breathvector,holconvector},'model','interaction','varnames',{'passive-active','hol-conn'});
computecohensd(nofpassive, nanmean(co2means(breathvector==1)), nanstd(co2means(breathvector==1)), nofactive, nanmean(co2means(breathvector==2)), nanstd(co2means(breathvector==2)))
computecohensd(nholo, nanmean(co2means(holconvector==1)), nanstd(co2means(holconvector==1)), nconcon, nanmean(co2means(holconvector==2)), nanstd(co2means(holconvector==2)))

[p,tbl,stats] = anovan(co2minimums,{breathvector,holconvector},'model','interaction','varnames',{'passive-active','hol-conn'});
computecohensd(nofpassive, nanmean(co2minimums(breathvector==1)), nanstd(co2minimums(breathvector==1)), nofactive, nanmean(co2minimums(breathvector==2)), nanstd(co2minimums(breathvector==2)))
computecohensd(nholo, nanmean(co2minimums(holconvector==1)), nanstd(co2minimums(holconvector==1)), nconcon, nanmean(co2minimums(holconvector==2)), nanstd(co2minimums(holconvector==2)))


%% Figures 1c and S1
clear p tbl stats ans

% Figure 1c
figure;
subplot (1,2,1)
hold on
errorbar(nanmean(co2times(breathvector==1,:)),nanstd(co2times(breathvector==1,:))/sqrt(ntotal),'Color', [0.5 0.5 0.5], 'LineWidth', 2)
errorbar(nanmean(co2times(breathvector==2,:)),nanstd(co2times(breathvector==2,:))/sqrt(ntotal),'k', 'LineWidth', 2)
xlim([0 7])
ylim([10 45])
hold off

subplot (1,2,2)
hold on
errorbar(nanmean(co2times(breathvector==1 & holconvector==1,:)),nanstd(co2times(breathvector==1& holconvector==1,:))/sqrt(nholo),'Color', [0.6 0.6 1], 'LineWidth', 2)
errorbar(nanmean(co2times(breathvector==2 & holconvector==1,:)),nanstd(co2times(breathvector==2& holconvector==1,:))/sqrt(nholo),'Color', [0.1 0.1 0.8], 'LineWidth', 2)
errorbar(nanmean(co2times(breathvector==1 & holconvector==2,:)),nanstd(co2times(breathvector==1& holconvector==2,:))/sqrt(nconcon),'Color', [1 0.6 0.6], 'LineWidth', 2)
errorbar(nanmean(co2times(breathvector==2 & holconvector==2,:)),nanstd(co2times(breathvector==2& holconvector==2,:))/sqrt(nconcon),'Color', [0.8 0.1 0.1], 'LineWidth', 2)
xlim([0 7])
ylim([10 45])
hold off


% Figure S1
figure;
subplot (1,2,1)
plot(co2times(breathvector==1 & holconvector==1,:)','Color', [0.6 0.6 1], 'LineWidth', 2)
xlim([0 7])
ylim([0 60])
subplot (1,2,2)
plot(co2times(breathvector==1 & holconvector==2,:)','Color', [1 0.6 0.6], 'LineWidth', 2)
xlim([0 7])
ylim([0 60])

figure;
subplot (1,2,1)
plot(co2times(breathvector==2 & holconvector==1,:)','Color', [0.1 0.1 0.8], 'LineWidth', 2)
xlim([0 7])
ylim([0 60])
subplot (1,2,2)
plot(co2times(breathvector==2 & holconvector==2,:)','Color', [0.8 0.1 0.1], 'LineWidth', 2)
xlim([0 7])
ylim([0 60])


%% Figure 2A and S2
%DASCnames = {'Unity'; 'Spiritual Experience';'Bliss';'Insight';'Disembodiment';'Impaired Control';'Anxiety'; 'Complex Imagery'; 'Simple Imagery'; 'Synaesthesia'; 'Changed Meaning' };

% Figure 2a
figure;
hold on

breathcolor = [1 0.5 1;0.7 0.2 0.7];
for br = 1:2
    DASCtemp = DASCmatrix(breathvector==br,:);
    plot(nanmean(DASCtemp),'Color',breathcolor(br,:), 'LineWidth', 2)         
end


xticks([1:11])
xticklabels({'Unity','Spirit','Bliss','Insight','Disembody','Impaired Control','Anxiety','Complex Image','Elementary Image','Synesthesia','Changed Meaning'})
xtickangle(90)
ylim([0 65])
xlim([0 12])


plot (psiloDASC,'--k','LineWidth',1)
plot (mdmaDASC,'--','Color', [0.5 0.5 0.5], 'LineWidth',1)
plot (lsdDASC,'-.k','LineWidth',1)
plot (placeboDASC,'-','Color', [0.5 0.5 0.5],'LineWidth',1)

hold off


% Extra figure - breathwork styles
figure;
hold on

techcolor = {[0.6 0.6 1;1 0.6 0.6],[0.1 0.1 0.8 ; 0.8 0.1 0.1]};
% Holo in blue, Conn in red, lighter colors = passive breath

for tech = 1:2
    for br = 1:2
        DASCtemp = DASCmatrix(breathvector==br&holconvector==tech,:) 
        errorbar(nanmean(DASCtemp),nanstd(DASCtemp)./sqrt(size(DASCtemp,1)),'Color',techcolor{br}(tech,:), 'LineWidth', 2)         
    end
end

xticks([1:11])
xticklabels({'Unity','Spirit','Bliss','Insight','Disembody','Impaired Control','Anxiety','Complex Image','Elementary Image','Synesthesia','Changed Meaning'})
xtickangle(90)
ylim([0 55])
xlim([0 12])
hold off


DASCtemp = DASCmatrix(breathvector==2,:); 
cohensdtemp = (nanmean(DASCtemp)-placeboDASC)./nanstd(DASCtemp)
cohensdtemp = (nanmean(DASCtemp)-psiloDASC)./nanstd(DASCtemp)
cohensdtemp = (nanmean(DASCtemp)-lsdDASC)./nanstd(DASCtemp)
cohensdtemp = (nanmean(DASCtemp)-mdmaDASC)./nanstd(DASCtemp)

DASCtemp = DASCmatrix(breathvector==1,:); 
cohensdtemp = (nanmean(DASCtemp)-placeboDASC)./nanstd(DASCtemp)



for br = 1:2
    DASCtemp = DASCmatrix(breathvector==br,:); 
    
    [h,p,ci,stats] = ttest(DASCtemp-psiloDASC);
    psiloDASCttestP{br} = p;
    psiloDASCttestStats{br} = stats;
    psiloDASCttestH{br} = p < (1-0.95^(1/11));% sidak correction family-wise error rate
    
    [h,p,ci,stats] = ttest(DASCtemp-mdmaDASC);
    mdmaDASCttestP{br} = p;
    mdmaDASCttestStats{br} = stats;
    mdmaDASCttestH{br} = p < (1-0.95^(1/11));% sidak correction family-wise error rate

    [h,p,ci,stats] = ttest(DASCtemp-lsdDASC);
    lsdDASCttestP{br} = p;
    lsdDASCttestStats{br} = stats;
    lsdDASCttestH{br} = p < (1-0.95^(1/11));% sidak correction family-wise error rate

    [h,p,ci,stats] = ttest(DASCtemp-placeboDASC);
    placeboDASCttestP{br} = p;
    placeboDASCttestStats{br} = stats;
    placeboDASCttestH{br} = p < (1-0.95^(1/11));
end



anovascaletemp = [];
anovabreathtemp = [];
anovaDASCtemp = [];

for scale = 1:11
    anovascaletemp = [anovascaletemp;ones(61,1)*scale];
    anovabreathtemp = [anovabreathtemp;breathvector];
    anovaDASCtemp = [anovaDASCtemp;DASCmatrix(:,scale)];
end

[p,tbl,stats] = anovan(anovaDASCtemp,{anovabreathtemp,anovascaletemp},'model','interaction','varnames',{'passive-active','DASC scale'});

% Cohen's D of main effects
vector1 = DASCmatrix(breathvector==1,:);
vector2 = DASCmatrix(breathvector==2,:);
vector1 = vector1(:);
vector2 = vector2(:);
computecohensd(length(vector1),nanmean(vector1),nanstd(vector1),length(vector2),nanmean(vector2),nanstd(vector2))

cohensdtemp = [];
for subs1 = 1:10
    for subs2 = subs1+1:11
        
        vector1 = DASCmatrix(:,subs1);
        vector2 = DASCmatrix(:,subs2);
        vector1 = vector1(:);
        vector2 = vector2(:);
        cohensdtemp(subs1,subs2) = computecohensd(length(vector1),nanmean(vector1),nanstd(vector1),length(vector2),nanmean(vector2),nanstd(vector2));
    end
end

hold off

clear h p ci tbl stats anovascaletemp anovabreathtemp anovaDASCtemp scale
clear vector1 vector2 subs1 subs2 cohensdtemp
clear tech br

%% Figure 2B
clear p tbl stats *ttest* ans
% techcolor = {[0.6 0.6 1;1 0.6 0.6],[0.1 0.1 0.8 ; 0.8 0.1 0.1]};
% Holo in blue, Conn in red, lighter colors = passive breath


% Figure 2b
figure;
hold on
breathcolor = [1 0.5 1;0.7 0.2 0.7];

for br = 1:2
    MEQtemp = MEQmatrix(breathSURVEYvector==br,1:4); 
    plot(nanmean(MEQtemp),'Color',breathcolor(br,:), 'LineWidth', 2)

end


xticks([1:4])
xticklabels({'Transcendent','Pos. Mood','Ineffable','Mystical'})
xtickangle(90)
ylim([0 1])
xlim([0 5])

plot (psiloMEQ,'--k','LineWidth',1)
plot (lsdMEQ,'-.k','LineWidth',1)
plot (mdmaMEQ,'--','Color', [0.5 0.5 0.5],'LineWidth',1)
plot (placeboMEQ,'Color', [0.5 0.5 0.5],'LineWidth',1)

hold off


MEQtemp = MEQmatrix(breathSURVEYvector==2,1:4) 
cohensdtemp = nanmean(MEQtemp-placeboMEQ)./nanstd(MEQtemp)
cohensdtemp = nanmean(MEQtemp-psiloMEQ)./nanstd(MEQtemp)
cohensdtemp = nanmean(MEQtemp-lsdMEQ)./nanstd(MEQtemp)
cohensdtemp = nanmean(MEQtemp-mdmaMEQ)./nanstd(MEQtemp)
clear MEQtemp cohensdtemp

%Same figure, but step by step
figure;
breathcolor = [1 0.5 1;0.7 0.2 0.7];

subplot(1,4,1)
hold on
plot (placeboMEQ,'Color', [0.5 0.5 0.5],'LineWidth',1)
xticks([1:4])
xticklabels({'Transcendent','Pos. Mood','Ineffable','Mystical'})
xtickangle(90)
ylim([0 1])
xlim([0 5])


subplot(1,4,2)
hold on
plot (placeboMEQ,'Color', [0.5 0.5 0.5],'LineWidth',1)
plot (psiloMEQ,'--k','LineWidth',1)
plot (lsdMEQ,'-.k','LineWidth',1)
plot (mdmaMEQ,'--','Color', [0.5 0.5 0.5],'LineWidth',1)
xticks([1:4])
xticklabels({'Transcendent','Pos. Mood','Ineffable','Mystical'})
xtickangle(90)
ylim([0 1])
xlim([0 5])


subplot(1,4,3)
hold on
plot (placeboMEQ,'Color', [0.5 0.5 0.5],'LineWidth',1)
plot (psiloMEQ,'--k','LineWidth',1)
plot (lsdMEQ,'-.k','LineWidth',1)
plot (mdmaMEQ,'--','Color', [0.5 0.5 0.5],'LineWidth',1)
MEQtemp = MEQmatrix(breathSURVEYvector==2,1:4); 
plot(nanmean(MEQtemp),'Color',breathcolor(2,:), 'LineWidth', 2)
xticks([1:4])
xticklabels({'Transcendent','Pos. Mood','Ineffable','Mystical'})
xtickangle(90)
ylim([0 1])
xlim([0 5])


subplot(1,4,4)
hold on
plot (placeboMEQ,'Color', [0.5 0.5 0.5],'LineWidth',1)
plot (psiloMEQ,'--k','LineWidth',1)
plot (lsdMEQ,'-.k','LineWidth',1)
plot (mdmaMEQ,'--','Color', [0.5 0.5 0.5],'LineWidth',1)
for br = 1:2
    MEQtemp = MEQmatrix(breathSURVEYvector==br,1:4); 
    plot(nanmean(MEQtemp),'Color',breathcolor(br,:), 'LineWidth', 2)

end
xticks([1:4])
xticklabels({'Transcendent','Pos. Mood','Ineffable','Mystical'})
xtickangle(90)
ylim([0 1])
xlim([0 5])




% Extra figure - breathwork styles

figure;
hold on

for tech = 1:2
    for br = 1:2
        MEQtemp = MEQmatrix(breathSURVEYvector==br&holconSURVEYvector==tech,1:4); 
        errorbar(nanmean(MEQtemp),nanstd(MEQtemp)./sqrt(size(MEQtemp,1)),'Color',techcolor{br}(tech,:), 'LineWidth', 2)  
    end
end

xticks([1:4])
xticklabels({'Transcendent','Pos. Mood','Ineffable','Mystical'})
xtickangle(90)
ylim([0 1])
xlim([0 5])
hold off



clear placeboMEQttest* psiloMEQttest*
for br = 1:2
    MEQtemp = MEQmatrix(breathSURVEYvector==br,1:4); 
    [h,p,ci,stats] = ttest(MEQtemp-placeboMEQ);
    placeboMEQttestP{br} = p;
    placeboMEQttestStats{br} = stats;
    placeboMEQttestH{br} = p < (1-0.95^(1/4)); 

    [h,p,ci,stats] = ttest(MEQtemp-psiloMEQ);
    psiloMEQttestP{br} = p;
    psiloMEQttestStats{br} = stats;
    psiloMEQttestH{br} = p < (1-0.95^(1/4)); 
    
    [h,p,ci,stats] = ttest(MEQtemp-lsdMEQ);
    lsdMEQttestP{br} = p;
    lsdMEQttestStats{br} = stats;
    lsdMEQttestH{br} = p < (1-0.95^(1/4)); 
    
    [h,p,ci,stats] = ttest(MEQtemp-mdmaMEQ);
    mdmaMEQttestP{br} = p;
    mdmaMEQttestStats{br} = stats;
    mdmaMEQttestH{br} = p < (1-0.95^(1/4)); 
    
end

MEQtemp = MEQmatrix(breathSURVEYvector==1,1:4); 
cohensdtemp = (nanmean(MEQtemp)-placeboMEQ)./nanstd(MEQtemp)


anovascaletemp = [];
anovabreathtemp = [];
anovaMEQtemp = [];
for scale = 1:4
    anovascaletemp = [anovascaletemp;ones(33,1)*scale];
    anovabreathtemp = [anovabreathtemp;breathSURVEYvector];
    anovaMEQtemp = [anovaMEQtemp;MEQmatrix(:,scale)];
end

[p,tbl,stats] = anovan(anovaMEQtemp,{anovabreathtemp,anovascaletemp},'model','interaction','varnames',{'passive-active','MEQ scale'});

% Cohen's D of main effects
vector1 = MEQmatrix(breathSURVEYvector==1,:);
vector2 = MEQmatrix(breathSURVEYvector==2,:);
vector1 = vector1(:);
vector2 = vector2(:);
computecohensd(length(vector1),nanmean(vector1),nanstd(vector1),length(vector2),nanmean(vector2),nanstd(vector2))

cohensdtemp = [];
for subs1 = 1:3
    for subs2 = subs1+1:4
        
        vector1 = MEQmatrix(:,subs1);
        vector2 = MEQmatrix(:,subs2);
        vector1 = vector1(:);
        vector2 = vector2(:);
        cohensdtemp(subs1,subs2) = computecohensd(length(vector1),nanmean(vector1),nanstd(vector1),length(vector2),nanmean(vector2),nanstd(vector2));
    end
end

clear br tech anovascaletemp anovabreathtemp *MEQtemp h ci ans
clear scale subs* vector1 vector2 

%% Figure 2C-D
clear p  tbl stats *ttest*
close all

depthmeans = nanmean(depthtimes(:,sessionpeakwindow)')';
depthmaximums = max(depthtimes(:,sessionpeakwindow)')';

% for activepassive = 1:2
%     nanmean(depthmeans(breathvector == activepassive))
%     nanstd(depthmeans(breathvector == activepassive))/sqrt(length(find(breathvector == activepassive)))
%     nanmean(depthmaximums(breathvector == activepassive))
%     nanstd(depthmaximums(breathvector == activepassive))/sqrt(length(find(breathvector == activepassive)))
% end
% 
% 
% for tech = 1:2
%     nanmean(depthmeans(holconvector == tech & breathvector==1))
%     nanstd(depthmeans(holconvector == tech& breathvector==1))/sqrt(length(find(holconvector == tech & breathvector==1)))
%     nanmean(depthmaximums(holconvector == tech & breathvector==1))
%     nanstd(depthmaximums(holconvector == tech& breathvector==1))/sqrt(length(find(holconvector == tech & breathvector==1)))
% end

% Figure 2c-d
figure;
subplot(1,2,1)
boxplot(depthmeans, [breathvector,holconvector])
ylim([0 5])
subplot(1,2,2)
boxplot(depthmaximums, [breathvector,holconvector])
ylim([0 5])

[p,tbl,stats] = anovan(depthmeans,{breathvector,holconvector},'model','interaction','varnames',{'passive-active','hol-conn'});
[p,tbl,stats] = anovan(depthmaximums,{breathvector,holconvector},'model','interaction','varnames',{'passive-active','hol-conn'});

% Cohen's D comparisons for anovas
computecohensd(nofpassive, nanmean(depthmeans(breathvector==1)), nanstd(depthmeans(breathvector==1)), nofactive, nanmean(depthmeans(breathvector==2)), nanstd(depthmeans(breathvector==2)))
computecohensd(nholo, nanmean(depthmeans(holconvector==1)), nanstd(depthmeans(holconvector==1)), nconcon, nanmean(depthmeans(holconvector==2)), nanstd(depthmeans(holconvector==2)))
computecohensd(nofpassive, nanmean(depthmaximums(breathvector==1)), nanstd(depthmaximums(breathvector==1)), nofactive, nanmean(depthmaximums(breathvector==2)), nanstd(depthmaximums(breathvector==2)))
computecohensd(nholo, nanmean(depthmaximums(holconvector==1)), nanstd(depthmaximums(holconvector==1)), nconcon, nanmean(depthmaximums(holconvector==2)), nanstd(depthmaximums(holconvector==2)))

computecohensd(nofpassive, nanmean(depthmeans(breathvector==1 & holconvector==1)), nanstd(depthmeans(breathvector==1& holconvector==1)), nofactive, nanmean(depthmeans(breathvector==1 & holconvector==2)), nanstd(depthmeans(breathvector==1 & holconvector==2)))

% Extra figure  - individual trajectories
figure;
subplot (1,2,1)
plot(depthtimes(breathvector==2 & holconvector==1,:)','Color', [0.8 0.2 0.2], 'LineWidth', 2)
xlim([0 7])
ylim([0 5])
subplot (1,2,2)
plot(depthtimes(breathvector==2 & holconvector==2,:)','Color', [0.2 0.2 0.8], 'LineWidth', 2)
xlim([0 7])
ylim([0 5])

clear ans p tbl stats

%% Figure 2E
close all

techcolor = {[0.6 0.6 1;1 0.6 0.6],[0.1 0.1 0.8 ; 0.8 0.1 0.1]};
% Holo in blue, Conn in red, lighter colors = passive breath

figure;
subplot (1,3,1)
hold on
errorbar(nanmean(depthtimes(breathvector==1,:)),nanstd(depthtimes(breathvector==1,:))/sqrt(ntotal),'Color', [0.5 0.5 0.5], 'LineWidth', 2)
errorbar(nanmean(depthtimes(breathvector==2,:)),nanstd(depthtimes(breathvector==2,:))/sqrt(ntotal),'k', 'LineWidth', 2)
xlim([0 7])
ylim([0 5])
hold off

subplot (1,3,2)
hold on
errorbar(nanmean(depthtimes(breathvector==1 & holconvector==1,:)),nanstd(depthtimes(breathvector==1& holconvector==1,:))/sqrt(nholo),'Color', [1 0.6 0.6], 'LineWidth', 2)
errorbar(nanmean(depthtimes(breathvector==2 & holconvector==1,:)),nanstd(depthtimes(breathvector==2& holconvector==1,:))/sqrt(nholo),'Color', [0.8 0.1 0.1], 'LineWidth', 2)
errorbar(nanmean(depthtimes(breathvector==1 & holconvector==2,:)),nanstd(depthtimes(breathvector==1& holconvector==2,:))/sqrt(nconcon),'Color', [0.6 0.6 1], 'LineWidth', 2)
errorbar(nanmean(depthtimes(breathvector==2 & holconvector==2,:)),nanstd(depthtimes(breathvector==2& holconvector==2,:))/sqrt(nconcon),'Color', [0.1 0.1 0.8], 'LineWidth', 2)
xlim([0 7])
ylim([0 5])
hold off


subplot (1,3,3)
hold on
errorbar(nanmean(depthtimes(62:end,:)),nanstd(depthtimes(62:end,:))/sqrt(30),'Color', [1 0.8 0.8], 'LineWidth', 2)
errorbar(nanmean(depthtimes(breathvector==1 & holconvector==1,:)),nanstd(depthtimes(breathvector==1& holconvector==1,:))/sqrt(30),'Color', [1 0.6 0.6], 'LineWidth', 2)
errorbar(nanmean(depthtimes(breathvector==2 & holconvector==1,:)),nanstd(depthtimes(breathvector==2& holconvector==1,:))/sqrt(30),'Color', [0.8 0.1 0.1], 'LineWidth', 2)
errorbar(nanmean(depthtimes(breathvector==1 & holconvector==2,:)),nanstd(depthtimes(breathvector==1& holconvector==2,:))/sqrt(31),'Color', [0.6 0.6 1], 'LineWidth', 2)
errorbar(nanmean(depthtimes(breathvector==2 & holconvector==2,:)),nanstd(depthtimes(breathvector==2& holconvector==2,:))/sqrt(31),'Color', [0.1 0.1 0.8], 'LineWidth', 2)
xlim([0 7])
ylim([0 5])
hold off



%% Supp. Figure S2

clear r p tbl tech
close all

% depthmeans = nanmean(depthtimes(:,sessionpeakwindow)')';
% depthmaximums = max(depthtimes(:,sessionpeakwindow)')';


% DASC correlation with hand signs
DASCdepthcorrels = [];

for tech = 1:2
    DASCholcon = DASCmatrix(holconvector==tech,:);
    depthholcon = depthmeans(holconvector==tech);
    [r{tech},p{tech}]=corrcoef([DASCholcon,depthholcon]);
    DASCdepthcorrels(:,tech) = r{tech}(:,end-1);
    DASCdepthcorrelpvalues(:,tech) = p{tech}(:,end-1);
end

[r{3},p{3}]=corrcoef([DASCmatrix,depthmeans]);
DASCdepthcorrels(:,3) = r{3}(:,end);
DASCdepthcorrelpvalues(:,3) = p{3}(:,end);
DASCdepthcorrels (12,:)=[];
DASCdepthcorrelpvalues (12,:)=[];



nanmean(DASCdepthcorrels(:,3))
nanstd(DASCdepthcorrels(:,3))



% MEQ correlations with hand signs

MEQhandsignvector = depthmeans([1:15,31:48]);

for scale = 1:5
    [rMEQ,pMEQ]=corrcoef(MEQmatrix(:,scale),MEQhandsignvector,'rows','pairwise')
    MEQcorr(scale) = rMEQ(2);
    MEQcorrP(scale) = pMEQ(2);    
end

% total figure
figure;
hold on
scatter (ones(11,1)+rand(11,1)*0.2-0.1, DASCdepthcorrels(:,3),'MarkerFaceColor','k','MarkerEdgeColor','k')
plot([0.8 1.2],[1 1]*nanmean(DASCdepthcorrels(:,3)),'k','LineWidth', 2)

scatter (1+ones(5,1)+rand(5,1)*0.2-0.1, MEQcorr,'MarkerFaceColor','k','MarkerEdgeColor','k')
plot([1.8 2.2],[1 1]*nanmean(MEQcorr),'k','LineWidth', 2)


hold off
xlim ([0.5 2.5])
ylim([0 1])


clear scale tech r p rMEQ pMEQ DASCdepth*
%% Figure 3B
% Figure 3A computed in R script

techcolor = {[0.6 0.6 1;1 0.6 0.6],[0.1 0.1 0.8 ; 0.8 0.1 0.1]};
% Holo in blue, Conn in red, lighter colors = passive breath

[corrtemp,corrtempP] = corrcoef([co2means,depthmeans])
[corrtempH,corrtempHP] = corrcoef([co2means(holconvector==1),depthmeans(holconvector==1)],'rows','pairwise')
[corrtempC,corrtempCP] = corrcoef([co2means(holconvector==2),depthmeans(holconvector==2)],'rows','pairwise')

figure;
hold on
for tech = 1:2
    for br = 1:2
        if br == 1
            scatter(co2means(breathvector == br & holconvector == tech ),depthmeans(breathvector == br & holconvector == tech),'MarkerEdgeColor',techcolor{br}(tech,:))
        else
            scatter(co2means(breathvector == br & holconvector == tech ),depthmeans(breathvector == br & holconvector == tech),'MarkerEdgeColor',techcolor{br}(tech,:),'MarkerFaceColor',techcolor{br}(tech,:))
        end
    end
end

text(1,1,['r (holo) = ',num2str(corrtempH(2))])
text(1,0.5,['r (conn) = ',num2str(corrtempC(2))])
ylim([0 5])
xlim([0 50])



% Extra figure: raw data correlations across all measurement time points

corrtemp = corrcoef([co2times(:),depthtimes(:)],'rows','pairwise')

ctemp1 = co2times(holconvector==1,:);
ctemp2 = depthtimes(holconvector==1,:);
[corrtempH,corrtempHP]  = corrcoef([ctemp1(:),ctemp2(:)],'rows','pairwise')

ctemp1 = co2times(holconvector==2,:);
ctemp2 = depthtimes(holconvector==2,:);
[corrtempC,corrtempCP] = corrcoef([ctemp1(:),ctemp2(:)],'rows','pairwise')

figure;
hold on

for tech = 1:2
    for br = 1:2
        stemp1 = co2times(breathvector == br & holconvector == tech , :);
        stemp2 = depthtimes(breathvector == br & holconvector == tech, :);
        scatter(stemp1(:),stemp2(:),'MarkerEdgeColor',techcolor{br}(tech,:))%,'MarkerFaceColor',techcolor{br}(tech,:))
    end
end

text(1,5.8,['r (holo) = ',num2str(corrtempH(2))])
text(1,5.3,['r (conn) = ',num2str(corrtempC(2))])
ylim([0 6])
xlim([0 60])
hold off

clear tech br corrtemp* stemp* ctemp* ans p r

%% Figure 3C-D

% Figure 3c - Holotropic Trajectories
cmap = [0 0 0.6; 0 0 0.8; 0.1 0.1 0.9; 0.2 0.2 1 ; 0.5 0.5 1;0.7 0.7 1];

for br = 1:3
    co2MeanTrajectoriesHOLO(br,1:noftimepoints)= nanmean(co2times(find(breathvector==br-1 & holconvector==1),:));
    depthMeanTrajectoriesHOLO(br,1:noftimepoints)= nanmean(depthtimes(find(breathvector==br-1 & holconvector==1),:));
    co2SEMTrajectoriesHOLO(br,1:noftimepoints)= nanstd(co2times(find(breathvector==br-1 & holconvector==1),:))/sqrt(length(find(breathvector==br-1 & holconvector==1)));
    depthSEMTrajectoriesHOLO(br,1:noftimepoints)= nanstd(depthtimes(find(breathvector==br-1 & holconvector==1),:))/sqrt(length(find(breathvector==br-1 & holconvector==1)));
end


figure;
subplot(1,2,1)
hold on
for time = 1:noftimepoints
    % plot means
    plot(co2MeanTrajectoriesHOLO(2,time),depthMeanTrajectoriesHOLO(2,time),'o','MarkerFaceColor', [1 1 1],'MarkerEdgeColor', cmap(time,:),'markersize',5)
    plot(co2MeanTrajectoriesHOLO(3,time),depthMeanTrajectoriesHOLO(3,time),'o','MarkerFaceColor', cmap(time,:),'MarkerEdgeColor', cmap(time,:),'markersize',7)  
    
    % plot errorbars
    plot([1 1]*co2MeanTrajectoriesHOLO(2,time),[depthMeanTrajectoriesHOLO(2,time)-depthSEMTrajectoriesHOLO(2,time),depthMeanTrajectoriesHOLO(2,time)+depthSEMTrajectoriesHOLO(2,time)],'Color', cmap(time,:),'LineWidth',1)
    plot([1 1]*co2MeanTrajectoriesHOLO(3,time),[depthMeanTrajectoriesHOLO(3,time)-depthSEMTrajectoriesHOLO(3,time),depthMeanTrajectoriesHOLO(3,time)+depthSEMTrajectoriesHOLO(3,time)],'Color', cmap(time,:),'LineWidth',1)

    plot([co2MeanTrajectoriesHOLO(2,time)-co2SEMTrajectoriesHOLO(2,time), co2MeanTrajectoriesHOLO(2,time)+co2SEMTrajectoriesHOLO(2,time)],[1 1]*depthMeanTrajectoriesHOLO(2,time),'Color', cmap(time,:),'LineWidth',1)
    plot([co2MeanTrajectoriesHOLO(3,time)-co2SEMTrajectoriesHOLO(3,time), co2MeanTrajectoriesHOLO(3,time)+co2SEMTrajectoriesHOLO(3,time)],[1 1]*depthMeanTrajectoriesHOLO(3,time),'Color', cmap(time,:),'LineWidth',1)

    %plot connecting lines
    if time < noftimepoints
       plot(co2MeanTrajectoriesHOLO(2,time:time+1),depthMeanTrajectoriesHOLO(2,time:time+1),'-','Color', cmap(time,:),'LineWidth',1)           
       plot(co2MeanTrajectoriesHOLO(3,time:time+1),depthMeanTrajectoriesHOLO(3,time:time+1),'-','Color', cmap(time,:),'LineWidth',2)           
    end
end
hold off
xlim([0 50])
ylim([0 5])


% Figure 3d - Conscious-Connected Trajectories
cmap = [0.6 0 0; 0.8 0 0; 0.9 0.1 0.1; 1 0.2 0.2 ;1 0.5 0.5;1 0.7 0.7];

for br = 1:2
    co2MeanTrajectoriesCONN(br,1:noftimepoints)= nanmean(co2times(find(breathvector==br & holconvector==2),:));
    depthMeanTrajectoriesCONN(br,1:noftimepoints)= nanmean(depthtimes(find(breathvector==br& holconvector==2),:));
    co2SEMTrajectoriesCONN(br,1:noftimepoints)= nanstd(co2times(find(breathvector==br & holconvector==2),:))/sqrt(length(find(breathvector==br & holconvector==2)));
    depthSEMTrajectoriesCONN(br,1:noftimepoints)= nanstd(depthtimes(find(breathvector==br & holconvector==2),:))/sqrt(length(find(breathvector==br & holconvector==2)));
end

subplot(1,2,2)
hold on
for time = 1:noftimepoints
   % plot means
    plot(co2MeanTrajectoriesCONN(1,time),depthMeanTrajectoriesCONN(1,time),'o','MarkerFaceColor', [1 1 1],'MarkerEdgeColor', cmap(time,:),'markersize',5)
    plot(co2MeanTrajectoriesCONN(2,time),depthMeanTrajectoriesCONN(2,time),'o','MarkerFaceColor', cmap(time,:),'MarkerEdgeColor', cmap(time,:),'markersize',7)
   
    % plot errorbars
    plot([1 1]*co2MeanTrajectoriesCONN(1,time),[depthMeanTrajectoriesCONN(1,time)-depthSEMTrajectoriesCONN(1,time),depthMeanTrajectoriesCONN(1,time)+depthSEMTrajectoriesCONN(1,time)],'Color', cmap(time,:),'LineWidth',1)
    plot([1 1]*co2MeanTrajectoriesCONN(2,time),[depthMeanTrajectoriesCONN(2,time)-depthSEMTrajectoriesCONN(2,time),depthMeanTrajectoriesCONN(2,time)+depthSEMTrajectoriesCONN(2,time)],'Color', cmap(time,:),'LineWidth',1)

    plot([co2MeanTrajectoriesCONN(1,time)-co2SEMTrajectoriesCONN(1,time), co2MeanTrajectoriesCONN(1,time)+co2SEMTrajectoriesCONN(1,time)],[1 1]*depthMeanTrajectoriesCONN(1,time),'Color', cmap(time,:),'LineWidth',1)
    plot([co2MeanTrajectoriesCONN(2,time)-co2SEMTrajectoriesCONN(2,time), co2MeanTrajectoriesCONN(2,time)+co2SEMTrajectoriesCONN(2,time)],[1 1]*depthMeanTrajectoriesCONN(2,time),'Color', cmap(time,:),'LineWidth',1)

    %plot connecting lines
    if time < noftimepoints
       plot(co2MeanTrajectoriesCONN(1,time:time+1),depthMeanTrajectoriesCONN(1,time:time+1),'-','Color', cmap(time,:),'LineWidth',1)           
       plot(co2MeanTrajectoriesCONN(2,time:time+1),depthMeanTrajectoriesCONN(2,time:time+1),'-','Color', cmap(time,:),'LineWidth',2)           
    end
    
end
hold off
xlim([0 50])
ylim([0 5])



% Extra figure - direct comparison of trajectories

figure;
hold on

% Holotropic trajectories
cmap = [0 0 0.6; 0 0 0.8; 0.1 0.1 0.9; 0.2 0.2 1 ; 0.5 0.5 1;0.7 0.7 1];
for br = 1:3
    co2MeanTrajectoriesHOLO(br,1:noftimepoints)= nanmean(co2times(find(breathvector==br-1 & holconvector==1),:));
    depthMeanTrajectoriesHOLO(br,1:noftimepoints)= nanmean(depthtimes(find(breathvector==br-1 & holconvector==1),:));
    co2SEMTrajectoriesHOLO(br,1:noftimepoints)= nanstd(co2times(find(breathvector==br-1 & holconvector==1),:))/sqrt(length(find(breathvector==br-1 & holconvector==1)));
    depthSEMTrajectoriesHOLO(br,1:noftimepoints)= nanstd(depthtimes(find(breathvector==br-1 & holconvector==1),:))/sqrt(length(find(breathvector==br-1 & holconvector==1)));
end

for time = 1:noftimepoints
    % plot means
    plot(co2MeanTrajectoriesHOLO(2,time),depthMeanTrajectoriesHOLO(2,time),'o','MarkerFaceColor', [1 1 1],'MarkerEdgeColor', cmap(time,:),'markersize',5)
    plot(co2MeanTrajectoriesHOLO(3,time),depthMeanTrajectoriesHOLO(3,time),'o','MarkerFaceColor', cmap(time,:),'MarkerEdgeColor', cmap(time,:),'markersize',7)  
    
    % plot errorbars
    plot([1 1]*co2MeanTrajectoriesHOLO(2,time),[depthMeanTrajectoriesHOLO(2,time)-depthSEMTrajectoriesHOLO(2,time),depthMeanTrajectoriesHOLO(2,time)+depthSEMTrajectoriesHOLO(2,time)],'Color', cmap(time,:),'LineWidth',1)
    plot([1 1]*co2MeanTrajectoriesHOLO(3,time),[depthMeanTrajectoriesHOLO(3,time)-depthSEMTrajectoriesHOLO(3,time),depthMeanTrajectoriesHOLO(3,time)+depthSEMTrajectoriesHOLO(3,time)],'Color', cmap(time,:),'LineWidth',1)

    plot([co2MeanTrajectoriesHOLO(2,time)-co2SEMTrajectoriesHOLO(2,time), co2MeanTrajectoriesHOLO(2,time)+co2SEMTrajectoriesHOLO(2,time)],[1 1]*depthMeanTrajectoriesHOLO(2,time),'Color', cmap(time,:),'LineWidth',1)
    plot([co2MeanTrajectoriesHOLO(3,time)-co2SEMTrajectoriesHOLO(3,time), co2MeanTrajectoriesHOLO(3,time)+co2SEMTrajectoriesHOLO(3,time)],[1 1]*depthMeanTrajectoriesHOLO(3,time),'Color', cmap(time,:),'LineWidth',1)

    %plot connecting lines
    if time < noftimepoints
       plot(co2MeanTrajectoriesHOLO(2,time:time+1),depthMeanTrajectoriesHOLO(2,time:time+1),'-','Color', cmap(time,:),'LineWidth',1)           
       plot(co2MeanTrajectoriesHOLO(3,time:time+1),depthMeanTrajectoriesHOLO(3,time:time+1),'-','Color', cmap(time,:),'LineWidth',2)           
    end
end


% Conscious-connected trajectories
cmap = [0.6 0 0; 0.8 0 0; 0.9 0.1 0.1; 1 0.2 0.2 ;1 0.5 0.5;1 0.7 0.7];

for br = 1:2
    co2MeanTrajectoriesCONN(br,1:noftimepoints)= nanmean(co2times(find(breathvector==br & holconvector==2),:));
    depthMeanTrajectoriesCONN(br,1:noftimepoints)= nanmean(depthtimes(find(breathvector==br& holconvector==2),:));
    co2SEMTrajectoriesCONN(br,1:noftimepoints)= nanstd(co2times(find(breathvector==br & holconvector==2),:))/sqrt(length(find(breathvector==br & holconvector==2)));
    depthSEMTrajectoriesCONN(br,1:noftimepoints)= nanstd(depthtimes(find(breathvector==br & holconvector==2),:))/sqrt(length(find(breathvector==br & holconvector==2)));
end


for time = 1:noftimepoints
   % plot means
    plot(co2MeanTrajectoriesCONN(1,time),depthMeanTrajectoriesCONN(1,time),'o','MarkerFaceColor', [1 1 1],'MarkerEdgeColor', cmap(time,:),'markersize',5)
    plot(co2MeanTrajectoriesCONN(2,time),depthMeanTrajectoriesCONN(2,time),'o','MarkerFaceColor', cmap(time,:),'MarkerEdgeColor', cmap(time,:),'markersize',7)
   
    % plot errorbars
    plot([1 1]*co2MeanTrajectoriesCONN(1,time),[depthMeanTrajectoriesCONN(1,time)-depthSEMTrajectoriesCONN(1,time),depthMeanTrajectoriesCONN(1,time)+depthSEMTrajectoriesCONN(1,time)],'Color', cmap(time,:),'LineWidth',1)
    plot([1 1]*co2MeanTrajectoriesCONN(2,time),[depthMeanTrajectoriesCONN(2,time)-depthSEMTrajectoriesCONN(2,time),depthMeanTrajectoriesCONN(2,time)+depthSEMTrajectoriesCONN(2,time)],'Color', cmap(time,:),'LineWidth',1)

    plot([co2MeanTrajectoriesCONN(1,time)-co2SEMTrajectoriesCONN(1,time), co2MeanTrajectoriesCONN(1,time)+co2SEMTrajectoriesCONN(1,time)],[1 1]*depthMeanTrajectoriesCONN(1,time),'Color', cmap(time,:),'LineWidth',1)
    plot([co2MeanTrajectoriesCONN(2,time)-co2SEMTrajectoriesCONN(2,time), co2MeanTrajectoriesCONN(2,time)+co2SEMTrajectoriesCONN(2,time)],[1 1]*depthMeanTrajectoriesCONN(2,time),'Color', cmap(time,:),'LineWidth',1)

    %plot connecting lines
    if time < noftimepoints
       plot(co2MeanTrajectoriesCONN(1,time:time+1),depthMeanTrajectoriesCONN(1,time:time+1),'-','Color', cmap(time,:),'LineWidth',1)           
       plot(co2MeanTrajectoriesCONN(2,time:time+1),depthMeanTrajectoriesCONN(2,time:time+1),'-','Color', cmap(time,:),'LineWidth',2)           
    end
    
end
hold off
xlim([0 50])
ylim([0 5])


clear cmap br time *Trajectories* 

%% Supp. Figure S1 - Individual trajectories
close all

subjectvector = [60:ntotal];

if holconvector(subjectvector) == 1
    cmap = [0 0 0.6; 0 0 0.8; 0.1 0.1 0.9; 0.2 0.2 1 ; 0.5 0.5 1;0.7 0.7 1];
else
    cmap = [0.6 0 0; 0.8 0 0; 0.9 0.1 0.1; 1 0.2 0.2 ;1 0.5 0.5;1 0.7 0.7];
end


% A:
for subj = subjectvector
    figure(subj);
    hold on
    
    if breathvector(subj)== 1
        for time = 1:noftimepoints
           % plot means
            plot(co2times(subj,time),depthtimes(subj,time),'o','MarkerFaceColor', [1 1 1],'MarkerEdgeColor', cmap(time,:),'markersize',7)

            %plot connecting lines
            if time < noftimepoints
               plot(co2times(subj,time:time+1),depthtimes(subj,time:time+1),'-','Color', cmap(time,:),'LineWidth',1.5)           
            end
        end
    else
      for time = 1:noftimepoints
           % plot means
            plot(co2times(subj,time),depthtimes(subj,time),'o','MarkerFaceColor', cmap(time,:),'MarkerEdgeColor', cmap(time,:),'markersize',7)

            %plot connecting lines
            if time < noftimepoints
               plot(co2times(subj,time:time+1),depthtimes(subj,time:time+1),'-','Color', cmap(time,:),'LineWidth',2.5)           
            end
        end        
        
    end
    
    hold off
    xlim([0 60])
    ylim([0 5])
end



% Extra figure: correlation of experience depth to minimal CO2 

techcolor = {[0.6 0.6 1;1 0.6 0.6],[0.1 0.1 0.8 ; 0.8 0.1 0.1]};
% Holo in blue, Conn in red, lighter colors = passive breath

corrtemp = corrcoef([co2minimums,depthmaximums])
[corrtempH,corrtempHP] = corrcoef([co2minimums(holconvector==1),depthmaximums(holconvector==1)],'rows','pairwise')
[corrtempC,corrtempCP] = corrcoef([co2minimums(holconvector==2),depthmaximums(holconvector==2)],'rows','pairwise')

figure;
hold on
for tech = 1:2
    for br = 1:2
        scatter(co2minimums(breathvector == br & holconvector == tech ),depthmaximums(breathvector == br & holconvector == tech),'MarkerEdgeColor',techcolor{br}(tech,:))%,'MarkerFaceColor',techcolor{br}(tech,:))
    end
end

text(1,1,['r (holo) = ',num2str(corrtempH(2))])
text(1,0.5,['r (conn) = ',num2str(corrtempC(2))])
ylim([0 5.5])
xlim([0 50])

clear subj tech time corrtemp* cohensdtemp *Trajectories* subjectvector

%% Figure 4

close all

techcolor = {[0.6 0.6 1;1 0.6 0.6],[0.1 0.1 0.8 ; 0.8 0.1 0.1]};


% Figure 4a-b: QIDS change
figure;
subplot(1,2,1)
hold on
plot (QIDSmatrix(breathSURVEYvector==2 & holconSURVEYvector==1,1:2)','Color',[0.1 0.1 0.8],'LineWidth',1)
plot (QIDSmatrix(breathSURVEYvector==2 & holconSURVEYvector==2,1:2)','Color',[0.8 0.1 0.1],'LineWidth',1)
errorbar(nanmean(QIDSmatrix(breathSURVEYvector==2,1:2)),nanstd(QIDSmatrix(breathSURVEYvector==2,1:2))/sqrt(length(find(breathSURVEYvector==2))),'k','LineWidth',2)
xlim([0.5 2.5])
ylim([0 15])

subplot(1,2,2)
hold on
plot (QIDSmatrix(breathSURVEYvector==1 & holconSURVEYvector==1,1:2)','Color',[0.6 0.6 1],'LineWidth',1)
plot (QIDSmatrix(breathSURVEYvector==1 & holconSURVEYvector==2,1:2)','Color',[1 0.6 0.6],'LineWidth',1)
errorbar(nanmean(QIDSmatrix(breathSURVEYvector==1,1:2)),nanstd(QIDSmatrix(breathSURVEYvector==1,1:2))/sqrt(length(find(breathSURVEYvector==1))),'Color',[0.5 0.5 0.5],'LineWidth',2)
xlim([0.5 2.5])
ylim([0 15])


% nanmean(QIDSmatrix(breathSURVEYvector==2,1:2))
% nanstd(QIDSmatrix(breathSURVEYvector==2,1:2))
% nanmean(QIDSmatrix(breathSURVEYvector==1,1:2))
% nanstd(QIDSmatrix(breathSURVEYvector==1,1:2))
[h,p,ci,stats] = ttest(QIDSmatrix(breathSURVEYvector==2,1),QIDSmatrix(breathSURVEYvector==2,2))
[h,p,ci,stats] = ttest(QIDSmatrix(breathSURVEYvector==1,1),QIDSmatrix(breathSURVEYvector==1,2))


% Figure 4c-d : WEMBWS changes
figure;
subplot(1,2,1)
hold on
plot (WEMBWSmatrix(breathSURVEYvector==2 & holconSURVEYvector==1,1:2)','Color',[0.1 0.1 0.8],'LineWidth',1)
plot (WEMBWSmatrix(breathSURVEYvector==2 & holconSURVEYvector==2,1:2)','Color',[0.8 0.1 0.1],'LineWidth',1)
errorbar(nanmean(WEMBWSmatrix(breathSURVEYvector==2,1:2)),nanstd(WEMBWSmatrix(breathSURVEYvector==2,1:2))/sqrt(length(find(breathSURVEYvector==2))),'k','LineWidth',2)
xlim([0.5 2.5])
ylim([30 70])

subplot(1,2,2)
hold on
plot (WEMBWSmatrix(breathSURVEYvector==1 & holconSURVEYvector==1,1:2)','Color',[0.6 0.6 1],'LineWidth',1)
plot (WEMBWSmatrix(breathSURVEYvector==1 & holconSURVEYvector==2,1:2)','Color',[1 0.6 0.6],'LineWidth',1)
errorbar(nanmean(WEMBWSmatrix(breathSURVEYvector==1,1:2)),nanstd(WEMBWSmatrix(breathSURVEYvector==1,1:2))/sqrt(length(find(breathSURVEYvector==1))),'Color',[0.5 0.5 0.5],'LineWidth',2)
xlim([0.5 2.5])
ylim([30 70])


% nanmean(WEMBWSmatrix(breathSURVEYvector==2,1:2))
% nanstd(WEMBWSmatrix(breathSURVEYvector==2,1:2))
% nanmean(WEMBWSmatrix(breathSURVEYvector==1,1:2))
% nanstd(WEMBWSmatrix(breathSURVEYvector==1,1:2))
[h,p,ci,stats] = ttest(WEMBWSmatrix(breathSURVEYvector==2,1),WEMBWSmatrix(breathSURVEYvector==2,2))
[h,p,ci,stats] = ttest(WEMBWSmatrix(breathSURVEYvector==1,1),WEMBWSmatrix(breathSURVEYvector==1,2))

clear br cmap

%% Extra figures for Figure 4E
% Figure 4E is generated in R
clear h p ci stats ans

%Predicted sustained changes
WEMBWSchange= WEMBWSmatrix(:,2)-WEMBWSmatrix(:,1);
QIDSchange= QIDSmatrix(:,2)-QIDSmatrix(:,1);

%Predictor variables
co2meansSURVEY = co2means([1:15,31:48]);
QIDS_DASC = DASCmatrix([1:15,31:48],:);
QIDS_DASC_OB = nanmean(QIDS_DASC(:,1:4)')';%oceanic boundlessness
QIDS_DASC_AED = nanmean(QIDS_DASC(:,5:7)')';%anxious ego dissolution
QIDS_DASC_VR = nanmean(QIDS_DASC(:,8:11)')';%visual restructuring
MEQtotalscore = MEQmatrix(:,5);

% QIDS
figure;
subplot(1,5,1)
scatter(co2meansSURVEY,QIDSchange,'k')
[r,p]= corrcoef([co2meansSURVEY QIDSchange],'rows','pairwise')
text (10,-13,['r = ',num2str(r(2))])
text (10,-15,['p = ',num2str(p(2))])

subplot(1,5,2)
scatter(MEQtotalscore,QIDSchange,'k')
[r,p]= corrcoef([MEQtotalscore QIDSchange],'rows','pairwise')
text (1,-13,['r = ',num2str(r(2))])
text (1,-15,['p = ',num2str(p(2))])

subplot(1,5,3)
scatter(QIDS_DASC_OB ,QIDSchange,'k')
[r,p]= corrcoef([QIDS_DASC_OB QIDSchange],'rows','pairwise')
text (1,-13,['r = ',num2str(r(2))])
text (1,-15,['p = ',num2str(p(2))])

subplot(1,5,4)
scatter(QIDS_DASC_AED,QIDSchange,'k')
[r,p]= corrcoef([QIDS_DASC_AED QIDSchange],'rows','pairwise')
text (1,-13,['r = ',num2str(r(2))])
text (1,-15,['p = ',num2str(p(2))])

subplot(1,5,5)
scatter(QIDS_DASC_VR ,QIDSchange,'k')
[r,p]= corrcoef([QIDS_DASC_VR QIDSchange],'rows','pairwise')
text (1,-13,['r = ',num2str(r(2))])
text (1,-15,['p = ',num2str(p(2))])


% WEMBWS
figure;
subplot(1,5,1)
scatter(co2meansSURVEY,WEMBWSchange,'k')
[r,p]= corrcoef([co2meansSURVEY WEMBWSchange],'rows','pairwise')
text (20,30,['r = ',num2str(r(2))])
text (20,27,['p = ',num2str(p(2))])

subplot(1,5,2)
scatter(MEQtotalscore,WEMBWSchange,'k')
[r,p]= corrcoef([MEQtotalscore WEMBWSchange],'rows','pairwise')
text (1,30,['r = ',num2str(r(2))])
text (1,27,['p = ',num2str(p(2))])


subplot(1,5,3)
scatter(QIDS_DASC_OB ,WEMBWSchange,'k')
[r,p]= corrcoef([QIDS_DASC_OB WEMBWSchange],'rows','pairwise')
text (1,30,['r = ',num2str(r(2))])
text (1,27,['p = ',num2str(p(2))])

subplot(1,5,4)
scatter(QIDS_DASC_AED,WEMBWSchange,'k')
[r,p]= corrcoef([QIDS_DASC_AED WEMBWSchange],'rows','pairwise')
text (1,30,['r = ',num2str(r(2))])
text (1,27,['p = ',num2str(p(2))])

subplot(1,5,5)
scatter(QIDS_DASC_VR ,WEMBWSchange,'k')
[r,p]= corrcoef([QIDS_DASC_VR WEMBWSchange],'rows','pairwise')
text (1,30,['r = ',num2str(r(2))])
text (1,27,['p = ',num2str(p(2))])

clear MEQtotalscore QIDS_DASC* r p br cmap

%% Figure 5

techcolor = {[0.6 0.6 1;1 0.6 0.6],[0.1 0.1 0.8 ; 0.8 0.1 0.1]};
% Holo in blue, Conn in red, lighter colors = passive breath

% Figure 5a-b: Amylase change
figure;
subplot(1,2,1)
hold on
plot (amylasematrix(breathMOLECULARvector==2 & holconMOLECULARvector==1,1:2)','Color',[0.1 0.1 0.8],'LineWidth',1)
plot (amylasematrix(breathMOLECULARvector==2 & holconMOLECULARvector==2,1:2)','Color',[0.8 0.1 0.1],'LineWidth',1)
errorbar(nanmean(amylasematrix(breathMOLECULARvector==2,1:2)),nanstd(amylasematrix(breathMOLECULARvector==2,1:2))/sqrt(length(find(breathMOLECULARvector==2))),'k','LineWidth',2)
xlim([0.5 2.5])
ylim([1 7])

subplot(1,2,2)
hold on
plot (amylasematrix(breathMOLECULARvector==1 & holconMOLECULARvector==1,1:2)','Color',[0.6 0.6 1],'LineWidth',1)
plot (amylasematrix(breathMOLECULARvector==1 & holconMOLECULARvector==2,1:2)','Color',[1 0.6 0.6],'LineWidth',1)
errorbar(nanmean(amylasematrix(breathMOLECULARvector==1,1:2)),nanstd(amylasematrix(breathMOLECULARvector==1,1:2))/sqrt(length(find(breathMOLECULARvector==1))),'Color',[0.5 0.5 0.5],'LineWidth',2)
xlim([0.5 2.5])
ylim([1 7])

[h,p,ci,stats] = ttest(amylasematrix(breathMOLECULARvector==2,1),amylasematrix(breathMOLECULARvector==2,2))


% Figure 5c-d: IL1b change
figure;
subplot(1,2,1)
hold on
plot (IL1bmatrix(breathMOLECULARvector==2 & holconMOLECULARvector==1,1:2)','Color',[0.1 0.1 0.8],'LineWidth',1)
plot (IL1bmatrix(breathMOLECULARvector==2 & holconMOLECULARvector==2,1:2)','Color',[0.8 0.1 0.1],'LineWidth',1)
errorbar(nanmean(IL1bmatrix(breathMOLECULARvector==2,1:2)),nanstd(IL1bmatrix(breathMOLECULARvector==2,1:2))/sqrt(length(find(breathMOLECULARvector==2))),'k','LineWidth',2)
xlim([0.5 2.5])
ylim([1 8])

subplot(1,2,2)
hold on
plot (IL1bmatrix(breathMOLECULARvector==1 & holconMOLECULARvector==1,1:2)','Color',[0.6 0.6 1],'LineWidth',1)
plot (IL1bmatrix(breathMOLECULARvector==1 & holconMOLECULARvector==2,1:2)','Color',[1 0.6 0.6],'LineWidth',1)
errorbar(nanmean(IL1bmatrix(breathMOLECULARvector==1,1:2)),nanstd(IL1bmatrix(breathMOLECULARvector==1,1:2))/sqrt(length(find(breathMOLECULARvector==1))),'Color',[0.5 0.5 0.5],'LineWidth',2)
xlim([0.5 2.5])
ylim([1 8])

[h,p,ci,stats] = ttest(IL1bmatrix(breathMOLECULARvector==2,1),IL1bmatrix(breathMOLECULARvector==2,2))


% Test differences between breathwork styles
AmylaseChange = ELISAtable.amylase_change;
IL1bChange = ELISAtable.IL1b_change;
[h,p,ci,stats] = ttest2(AmylaseChange(breathMOLECULARvector== 2 & holconMOLECULARvector== 1),AmylaseChange(breathMOLECULARvector== 2 & holconMOLECULARvector== 2 ))
[h,p,ci,stats] = ttest2(IL1bChange(breathMOLECULARvector== 2 & holconMOLECULARvector== 1),IL1bChange(breathMOLECULARvector== 2 & holconMOLECULARvector== 2 ))


clear h p ci stats


%% Methods: Subject Info Analysis
ReadinSubjectInfo

Age = SubjectInfo.age;
Gender = SubjectInfo.gender;
Workstatus = SubjectInfo.work_status;

% nanmean(Age)
% nanmean(Age(breathvector==1))
% nanmean(Age(breathvector==2))


%%
clear *temp*
