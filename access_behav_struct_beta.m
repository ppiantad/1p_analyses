%%
% to access BehaveData from inside behav_struct for trial filtering
titles = {'Trial','Block','ForceFree','bigSmall','stTime','choiceTime'...
    'collectionTime','shock','omission','omissionALL','WL','WSLScode','win_stay','lose_shift','lose_omit','smallRew','bigRew'};

iter = 0;
% These .mat files contain 4 files processed using AllAnimalFPanalysis
% with the Analysis_BL2initiation_fn_beta.m function. The files are as
% follows:

    %behav_struct
    %corrTable
    %ts1_struct
    %zall_struct

% Uncomment the particular file below to load the .mat file (will create a list of inputs later in a function to tidy this up)  
% Can uncomment and load diff. sessions to graph data from separate groups,
% etc (e.g. vmPFC RM Early REW 1.2 vs vmPFC RM Late REW 1.2)
    
% load('vmPFC_RM_EARLY_SINGLE_SESSION_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5.mat')

% load('vmPFC_RM_LATE_2_SESSIONS_ALL_EPOC__CHOICE_TRANGE_-33_43_BL_-10_-5.mat')

% load('BLA_RM_EARLY_SINGLE_SESSION_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5.mat')

% load('BLA_RM_LATE_2_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5.mat')

% load('vmPFC_RewMag_ALL_EPOC_Choice.mat')

% load('BLA_RewMag_ALL_EPOC_Choice.mat')

% load('vmPFC_RDT_ALL_EPOC_Choice.mat')

% load('BLA_RDT_ALL_EPOC_Choice.mat')

% load('vmPFC_RDT_ALL_EPOC__BL_-10_-5_Choice.mat')

% load('BLA_RDT_ALL_EPOC__BL_-10_-5_Choice.mat')

% load('D1_RewMag_ALL_TRANGE_-33_43_EPOC_Choice_BL_-3_0_BEFORE_INITIATION')

% load('D2_RewMag_ALL_TRANGE_-33_43_EPOC_Choice_BL_-3_0_BEFORE_INITIATION')

% load('D1_RM_EARLY_1_or_2_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5.mat')

load('D2_RM_EARLY_1_or_2_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5.mat')

% load('D1_RM_LATE_2_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5.mat')

% load('D2_RM_LATE_2_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5.mat')

% load('D1_RDT_ALL_TRANGE_-33_43_EPOC_Choice_BL_-3_0_BEFORE_INITIATION')

% load('D2_RDT_ALL_TRANGE_-33_43_EPOC_Choice_BL_-3_0_BEFORE_INITIATION')

% load('D1-OP-1_RDT_ALL_TRANGE_-33_43_EPOC_Choice_BL_-3_0_BEFORE_INITIATION')

% load('D1_RDT_2_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5.mat')

% load ('D2_RDT_2_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5.mat')
%%
clear zall_filtered; clear trials; clear ii; clear jj; clear ZallMean;

collect_times = [];
choice_times = [];
row = 1; 
for ii = 1:size(behav_struct,2)
%     arg1 = {'REW' 1.2}
%     arg2 = {'BLOCK' 3}
    BehavData = cell2mat(behav_struct(ii).data);
    BehavData = array2table(BehavData);
    BehavData.Properties.VariableNames = titles;
    
%     collect_times_all = [collect_times_all; {collect_times}];
    %Trim size of BehavData to match zall_struct (sometimes BehavData has a
    %few extra fake trials)
    keep = size(zall_struct.data{1, (ii)}, 1);
    BehavData = BehavData(1:keep, :);
    for kk = 1:size(BehavData, 1)
        if BehavData.stTime(kk) == BehavData.choiceTime(kk);
            BehavData(kk,:)=[];
        end
    end
    [BehavData, trials, varargin] = TrialFilter(BehavData,'REW',0.3);
    behav_struct_filtered(ii).data = BehavData;
    
     
    if isempty(BehavData)
        disp(['Animal has none of the specified trials'])
        ZallMean(ii,:) = NaN(1, size(zall_struct.data{1, (ii)},2));
        zerror(ii,:) = NaN(1, size(zall_struct.data{1, (ii)},2));
        choice_lat(ii) = NaN;
        collect_lat(ii) = NaN;

    continue
    else
    for jj = 1:size(trials,1)
        zall_filtered(ii).data(jj,:)= zall_struct.data{1, (ii)}(trials{jj},:);
        jj = jj+1;
    end
    jj = 1;
    if size(zall_filtered(ii).data, 1) == 1
%         disp(['Animal only has one trial, ignoring']) % comment out if you want to omit data if animal only has 1 trial of a given type
%         ZallMean(ii,:) = NaN(1, size(zall_struct.data{1, (ii)},2)); % comment out if you want to omit data if animal only has 1 trial of a given type
%         zerror(ii,:) = NaN(1, size(zall_struct.data{1, (ii)},2)); % comment out if you want to omit data if animal only has 1 trial of a given type
        ZallMean(ii,:) = zall_filtered(ii).data(jj,:);
        zerror(ii,:) = zeros(1, size(zall_struct.data{1, (ii)},2));
        
    elseif size(zall_filtered(ii).data,1) ~= 1
        ZallMean(ii,:) = mean(zall_filtered(ii).data);
        zerror(ii,:) = std(zall_filtered(ii).data)/size(zall_filtered(ii).data,1);
        choice_lat(ii) = mean(behav_struct_filtered([ii]).data.choiceTime - behav_struct_filtered([ii]).data.stTime)';
        collect_lat(ii) = mean(behav_struct_filtered([ii]).data.collectionTime - behav_struct_filtered([ii]).data.choiceTime)';
    end
%     clear zall_filtered;
    ts1 = cell2mat(ts1_struct.data);
%     collect_times(ii,:) = (behav_struct_filtered([ii]).data.collectionTime) - (behav_struct_filtered([ii]).data.stTime);
%     choice_times(ii,:) = (behav_struct_filtered([ii]).data.choiceTime) - (behav_struct_filtered([ii]).data.stTime);
     
    ii = ii+1;
    jj = 1;

    end
end    

[~,cols]=size(ZallMean);

for ii=1:cols
    zAvg(ii)=nanmean(ZallMean(:,ii));
    values=ZallMean(:,ii);
    stdDev(ii)=nanstd(values);
end


SEM=stdDev./(sqrt(size(ZallMean, 1)));
lo=zAvg-SEM;
hi=zAvg+SEM;

ind = ts1(1,:) > -10 & ts1(1,:) < 10;
ts1forprism = ts1(:,ind);
zAvgforprism = nanmean(ZallMean(:,ind));
semforprism = SEM(:,ind);


iter = iter+1;

AUC=[]; % cue, shock
for qq = 1:size(ZallMean,1);
    AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < -0 & ts1(1,:) > -2));
    AUC(qq,2)=trapz(ZallMean(qq,ts1(1,:) > 0 & ts1(1,:) < 2));
    qq=qq+1;
end
AUC_mean(:,iter) = nanmean(AUC);
AUC_std(:,iter) = nanstd(AUC)/sqrt(numel(AUC(:,1)));

%%
hold on;
figure(1)
hold on



% Plot vertical line at epoch onset, time = 0

%for loop for 5 tone retrieval- retrieval is 50 tones though??

% a=0;
% for ii=1:50
%     p1=patch([a a+30 a+30 a], [-5 -5 20 20], [.8 1 1], 'EdgeColor','none');
%     a=a+35;
% end 

% vline = plot([0 0],[min(lo) max(hi)],'--','LineWidth',2)

% y = ylim; % current y-axis limits
vline = plot([0 0],[min(lo) max(hi)]);
vline.Color = 'k';
vline.LineStyle = '--';
vline.LineWidth = 2;

xxx=[ts1, ts1(end:-1:1)];
yyy=[lo, hi(end:-1:1)];
hp= fill(xxx,yyy,[ .8 .8 .8]); hold on;
set(hp,'FaceColor', [ .8 .8 .8],'EdgeColor','none','FaceAlpha',0.8,'EdgeAlpha',0.8);hold on;
p3=plot(ts1, zAvg, 'k', 'LineWidth', 2,'DisplayName','z-score'); hold on;



% Make a legend
% legend([p3 hp],{'Block 1','SEM'}, 'AutoUpdate', 'off');

set(gca,'Layer','top')

%%
c = distinguishable_colors(numel(corrTable.animalNames));

figure(2);
set(gca, 'ColorOrder', c, 'NextPlot', 'replacechildren');

plot(ts1, ZallMean, 'LineWidth',2);

ID_legend = corrTable.animalNames;
legend(ID_legend);

%% Calculate AUC for different time windows, store in AUC array
% AUC=[]; % cue, shock
% for qq = 1:size(ZallMean,1);
%     AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < -0 & ts1(1,:) > -2));
%     AUC(qq,2)=trapz(ZallMean(qq,ts1(1,:) > 0 & ts1(1,:) < 5));
%     qq=qq+1;
% end
% AUC_mean(:,iter) = nanmean(AUC);
% AUC_std(:,iter) = nanstd(AUC)/sqrt(numel(AUC(:,1)));
%%
means=[]; % cue, shock
for qq = 1:size(ZallMean,1);
    means(qq,1)=mean(ZallMean(qq,ts1(1,:) < -0 & ts1(1,:) > -2));
    means(qq,2)=mean(ZallMean(qq,ts1(1,:) > 0 & ts1(1,:) < 5));
    qq=qq+1;
end
means_mean = mean(means);
means_std = std(means)/sqrt(numel(means(:,1)));

%%
%AUC correlation w/ choice latency
figure;
scatter(choice_lat,AUC(:,1));
coefficients = polyfit(choice_lat, AUC(:,1), 1);
xFit = linspace(0, max(choice_lat), 1000); %min(WSLS_Table.RiskPercent)
yFit = polyval(coefficients , xFit);
hold on;
title('Scatter: Choice Latency % vs Pre-Choice -2 to 0 AUC')
xlabel('Choice Latency')
ylabel('Normalized Pre-Choice AUC')
plot(xFit, yFit, 'r-', 'LineWidth', 2);
grid on;
choice_lat_AUC_corr = corrcoef(choice_lat, AUC(:,1),'rows','complete'); %correlate choice latency w/ pre-choice AUC
%%
%if a given animal has 0 trials, will need to delete rows to deal with
%mismatch
figure;
scatter(corrTable.RiskPercent,AUC(:,1));
coefficients = polyfit(corrTable.RiskPercent, AUC(:,1), 1);
xFit = linspace(0, max(corrTable.RiskPercent), 1000); %min(WSLS_Table.RiskPercent)
yFit = polyval(coefficients , xFit);
hold on;
title('Scatter: Risk % vs Pre-Choice')
xlabel('Proportion Risky Choices')
ylabel('Normalized Pre-Choice AUC')
plot(xFit, yFit, 'r-', 'LineWidth', 2);
grid on;

figure;
scatter(corrTable.RiskPercent,AUC(:,2));
coefficients = polyfit(corrTable.RiskPercent, AUC(:,2), 1);
xFit = linspace(0, max(corrTable.RiskPercent), 1000); %min(WSLS_Table.RiskPercent)
yFit = polyval(coefficients , xFit);
hold on;
title('Scatter: Risk % vs Post-Choice')
xlabel('Proportion Risky Choices')
ylabel('Normalized Post-Choice AUC')
plot(xFit, yFit, 'r-', 'LineWidth', 2);
grid on;

figure;
scatter(corrTable.LoseShiftPercent,AUC(:,1));
coefficients = polyfit(corrTable.LoseShiftPercent, AUC(:,1), 1);
xFit = linspace(0, max(corrTable.LoseShiftPercent), 1000); %min(WSLS_Table.RiskPercent)
yFit = polyval(coefficients , xFit);
hold on;
title('Scatter: LoseShift vs Pre-Choice')
xlabel('Proportion of Lose Shifts')
ylabel('Normalized Pre-Choice AUC')
plot(xFit, yFit, 'r-', 'LineWidth', 2);
grid on;

%%
figure;
hold on
C = [.8 .2 .8;
     .2 .75 .2];
CE = [.5 .1 .5];
superbar(AUC_mean(1,:), 'BarFaceColor', C, 'BarWidth', .5, 'E', AUC_std(1,:), 'ErrorbarColor', CE);
xlim([0.5 3.5]);

% title('superbar');

%% IF USING ALL SESSIONS + Need to double check that any mice with 1 session only are removed

LoseShift_D(:,1) = corrTable.LoseShiftPercent(1:2:end);
LoseShift_D(:,2) = corrTable.LoseShiftPercent(2:2:end);
WinStay_D(:,1) = corrTable.WinStayPercent(1:2:end);
WinStay_D(:,2) = corrTable.WinStayPercent(2:2:end);
RiskyPercent_D(:,1) = corrTable.RiskPercent(1:2:end);
RiskyPercent_D(:,2) = corrTable.RiskPercent(2:2:end);
meanRiskyPercent_D = mean(RiskyPercent_D);
SEMRiskyPercent_D = std(RiskyPercent_D)/sqrt(numel(RiskyPercent_D(:,1)));

%%
figure;
C = [.8 .2 .8;
     .2 .75 .2];
CE = [.5 .1 .5];
superbar(meanRiskyPercent_D, 'BarFaceColor', C, 'BarWidth', .5, 'E', SEMRiskyPercent_D, 'ErrorbarColor', CE);
xlim([0.5 2.5]);
% title('superbar');

%%
num_animals = size(zAvgforprism,1);

rc = num_animals/5;

sp_rows = floor(rc);
sp_cols = ceil(rc);

for i = 1:num_animals
    axi = subplot(sp_rows,sp_cols,i);
    plot(ts1forprism,zAvgforprism(i,:));
    title(ID_legend{i});
end

%%

for qq = 1:size(zall_filtered, 1)
    for pp = 1:size(zall_filtered(qq).data, 1)
        BL_shifted(pp,qq)=[behav_struct_filtered(qq).data.stTime(pp) - behav_struct_filtered(qq).data.choiceTime(pp)]; 
        peaks(pp, qq) = findchangepts(zall_filtered(qq).data(pp,ts1(1,:) > BL_shifted(pp) & ts1(1,:) < 0));    
        peaks_norm(pp, qq) = peaks(pp, qq) / abs(BL_shifted(pp,qq));
    end
    clear pp
end
%     ind = ts1(1,:) < BL_shifted(2) & ts1(1,:) > BL_shifted(1);
%     findchangepts(zall_filtered(qq).data(:,ts1(1,:) > 0 & ts1(1,:) < 2))
%%
for qq = 1:size(ZallMean,1)
    premean(qq,1)=mean(ZallMean(qq,ts1(1,:) < -0 & ts1(1,:) > -2));
    test_change(qq,1)=findchangepts(ZallMean(qq,ts1(1,:) < 0 & ts1(1,:) > -2));
    qq=qq+1;
end
increasing_means = ZallMean([1:6,13:14,19:20],:);
decreasing_means = ZallMean([7:12,15:18,21:26],:);
