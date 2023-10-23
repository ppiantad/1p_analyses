%% Changelog:
%   3/7/2019: Changed the way STREAM_STORES are read into MATLAB, because
%   different sessions have different stream names occasionally
%   3/18/2019: Works with Photometry System 2; Sensor C and D (as long as
%   they are named correctly)
%%
clear all; clc; close all;
%%
root = upper('h:');

BLOCKPATH =  strcat(root,'\MATLAB\TDTbin2mat\Photometry\RRD258\RRD258-220921-122157');

data = TDTbin2mat(BLOCKPATH);
%%
Sensor = input('Please enter name of sensor used (A,B, C, or D): \n','s'); %input which sensor the mouse was recorded from 
Sensor = upper(Sensor);

%find time of first TTL.  This value will be passed to Behavioral Data so
%that the timestamps of the ABET data will match the timestamps of Synapse.
%This does not work perfect for mice < #RRD47, because Sensor D was in Box
%1 (with TTL_), and Sensor C was in box 2 (with TL2_)
if Sensor == 'A';
    timeStart=data.epocs.TTL_.onset(1);
elseif Sensor == 'B';
    timeStart=data.epocs.TTL_.onset(1);  
elseif Sensor == 'C';
    timeStart=data.epocs.TTL_.onset(1);
elseif Sensor == 'D';
    timeStart=data.epocs.TL2_.onset(1);
elseif Sensor == 'RCAMP';
    timeStart=data.epocs.TTL_.onset(1);
end

%work behavioral data into photo data, generate descriptive statistics
[BehavData,ABETfile,Descriptives, block_end]=ABET2TableFn_Chamber_Av4('RRD179 04092021.csv',[]);


%align timestamps of behavioral data to timestamps of photometry by adding
%the time elapsed from when Synapse recording began to when the ABET
%program was started back to each relevant time column (choiceTime,
%collectionTime, stTime)
timeShift=timeStart*ones(numel(BehavData.choiceTime(:)),1);
BehavData.choiceTime(:)=BehavData.choiceTime(:)+timeShift;
BehavData.collectionTime(:)=BehavData.collectionTime(:)+timeShift;
BehavData.stTime(:)=BehavData.stTime(:)+timeShift;
block_end = block_end +timeShift;

%% added to attempt to integrate Approach/Aborts from BORIS

% **Add BORIS file here (e.g. '97 RDT 20200116_checked.csv', if it exists, otherwise leave empty**
boris_file = [];


if ~isempty(boris_file)
    [~,~,boris] = xlsread(boris_file);
    boris_extract = boris(17:end,[1,6]);
    boris_Extract_tbl = cell2table(boris_extract)
    
    boris_Extract_tbl.Properties.VariableNames = {'choiceTime','type'};
    
    boris_Extract_tbl.stTime = boris_Extract_tbl.choiceTime;
    boris_Extract_tbl.collectionTime = boris_Extract_tbl.choiceTime+30;
    
    % outerjoin(boris_Extract_tbl, BehavData);
    [boris_rows,~]=size(boris_Extract_tbl);
    for hh = 1:boris_rows
        if boris_Extract_tbl.stTime(hh) < block_end(1,1)
            boris_Extract_tbl.Block(hh) = 1;
        elseif boris_Extract_tbl.choiceTime(hh) > block_end(1,1) &&  boris_Extract_tbl.choiceTime(hh) < block_end(1,2)
            boris_Extract_tbl.Block(hh) = 2;
        elseif boris_Extract_tbl.choiceTime(hh) > block_end(1,2)
            boris_Extract_tbl.Block(hh,1) = 3;
        end
        
    end
    boris_Extract_tbl.Properties.VariableNames = {'choiceTime','type','stTime','collectionTime','Block'};
    abet_and_boris = tblvertcat(boris_Extract_tbl, BehavData);
    tblB = sortrows(abet_and_boris,{'choiceTime'},{'ascend'});
    
    BehavData = tblB;
end

%%
%pick trial types

BehavData=TrialFilter(BehavData,'REW',1.2);



%%
%create epoc for trials aligned to choice
data.epocs.Choice.name = 'Choice';
data.epocs.Choice.onset = BehavData.choiceTime; %
data.epocs.Choice.offset = BehavData.choiceTime+0.1*ones(numel(BehavData.choiceTime),1);
data.epocs.Choice.data = [1:numel(BehavData.choiceTime)]';


%create epoc for trials aligned to collection time
data.epocs.Collect.name = 'Collect';
data.epocs.Collect.onset = BehavData.collectionTime; %
data.epocs.Collect.offset = BehavData.collectionTime+0.1*ones(numel(BehavData.collectionTime),1);
data.epocs.Collect.data = [1:numel(BehavData.collectionTime)]';


%create epoc for trials aligned to Trial Start time
data.epocs.start.name = 'start';
data.epocs.start.onset = BehavData.stTime; %
data.epocs.start.offset = BehavData.stTime+0.1*ones(numel(BehavData.stTime),1);
data.epocs.start.data = [1:numel(BehavData.stTime)];

%%
%use TDT filter to select just the trials you want
EPOC = 'Choice'; %set which epoc to look at
streamNames = fieldnames(data.streams); 
% check for known names that we have used for 405 / 470 channels here,
% selects the one that is used for this particular session. if name
% stream names change, must be added to this loop manually
if Sensor == 'A';
    STREAM_STORE1 = 'x405S'
    STREAM_STORE2 = 'x470S'
elseif Sensor =='B';
    STREAM_STORE1 = 'x405G'
    STREAM_STORE2 = 'x470G'
elseif Sensor == 'D';
    STREAM_STORE1 = streamNames{2};%'d5sP' 
    STREAM_STORE2 = streamNames{1}; %'d0sP'
elseif Sensor == 'C';
    STREAM_STORE1 = streamNames{3};%'x45sP'
    STREAM_STORE2 = streamNames{4};%'x40sP'
elseif Sensor == 'RCAMP'
    STREAM_STORE1 = streamNames{3}; % 'x405C'
    STREAM_STORE2 = streamNames{4}; % 'x565G'
end
% STREAM_STORE1 = input('Please enter name of 405 stream: \n','s'); %check data.streams 
% STREAM_STORE2 = input('Please enter name of 470 stream: \n','s'); %check data.streams
% if ismember('x405G',StreamFields);
%     STREAM_STORE1 = 'x405G';
% elseif ismember('x405S',StreamFields);
%     STREAM_STORE1 = 'x405S';
% elseif ismember('x45sP',StreamFields);
%     STREAM_STORE1 = 'x45sP';
% end;
% 
% if ismember('x470G',StreamFields);
%     STREAM_STORE2 = 'x470G';
% elseif ismember('x470S',StreamFields);
%     STREAM_STORE2 = 'x470S';
% elseif ismember('x40sP',StreamFields);
%     STREAM_STORE2 = 'x40sP';
% end;
% if ~any(strcmp(StreamFields,'x405G'));
%      STREAM_STORE1 = 'x405G';
% elseif ~any(strcmp(StreamFields, 'x405S'));
%     STREAM_STORE1 = 'x405S';
% elseif ~any(strcmp(StreamFields, 'x45sP'));
%     STREAM_STORE1 = 'x45sP';
% end
% STREAM_STORE1 = isfields(data.streams('x405S')); %UV stream? can we not add to code to check if x405G / x405S / x45sP exist, and if so select that as the STREAM_STORE1 name?
% STREAM_STORE2 = isfields(data.streams('x470G','x470S','x45sP')); %GCamp stream?
%TRANGE sets the window size, data will be filtered around EPOC (whatever event), size of data will be window duration (2nd # in TRANGE) x sampling
%rate (1017) approximately
TRANGE = [-33 43]; % window size [start time relative to epoc onset, window duration] -10 20 default -4 12 -10 18 -31 41
BASELINE_PER = [-10 -5]; % baseline period within our window -6 -1 default; THE FIRST VALUE CAN'T BE < THE FIRST TRANGE VALUE -8 0 is ITI
ARTIFACT = Inf; % optionally set an artifact rejection level
data = TDTfilter(data,EPOC,'TIME',TRANGE); %perform TDTfilter 
%%
% if EPOC == 'Trial';
%     BL_time = BehavData.stTime(:)-BehavData.choiceTime(:);
%     time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);
% elseif EPOC == 'Collect';
%     BL_time = BehavData.collectionTime(:)-BehavData.choiceTime(:);
% end

time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);

if strcmpi('start',EPOC)
    time2EPOC = BehavData.stTime(:) - BehavData.choiceTime(:);

elseif strcmpi('Collect', EPOC)
    BL_time = BehavData.collectionTime(:)-BehavData.choiceTime(:);
    time2EPOC = BehavData.choiceTime(:)-BehavData.collectionTime(:); %time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);
elseif strcmpi('Choice', EPOC)
    BL_time = BehavData.stTime(:)-BehavData.choiceTime(:);
%     time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);
    time2EPOC = BehavData.stTime(:)-BehavData.choiceTime(:);
end

%%
% Optionally remove artifacts. If any waveform is above ARTIFACT level, or
% below -ARTIFACT level, remove it from the data set.
art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
good = ~art1 & ~art2;
data.streams.(STREAM_STORE1).filtered = data.streams.(STREAM_STORE1).filtered(good);

art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
good2 = ~art1 & ~art2;
data.streams.(STREAM_STORE2).filtered = data.streams.(STREAM_STORE2).filtered(good2);

numArtifacts = sum(~good) + sum(~good2);

%%
% Applying a time filter to a uniformly sampled signal means that the
% length of each segment could vary by one sample.  Let's find the minimum
% length so we can trim the excess off before calculating the mean.
minLength1 = min(cellfun('prodofsize', data.streams.(STREAM_STORE1).filtered));
minLength2 = min(cellfun('prodofsize', data.streams.(STREAM_STORE2).filtered));
data.streams.(STREAM_STORE1).filtered = cellfun(@(x) x(1:minLength1), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false);
data.streams.(STREAM_STORE2).filtered = cellfun(@(x) x(1:minLength2), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false);

allSignals = cell2mat(data.streams.(STREAM_STORE1).filtered');

% downsample 10x and average 405 signal
N = 10;
F405 = zeros(size(allSignals(:,1:N:end-N+1)));
for ii = 1:size(allSignals,1)
    F405(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
end
minLength1 = size(F405,2);

% Create mean signal, standard error of signal, and DC offset of 405 signal
meanSignal1 = mean(F405);
stdSignal1 = std(double(F405))/sqrt(size(F405,1));
dcSignal1 = mean(meanSignal1);

% downsample 10x and average 465 signal
allSignals = cell2mat(data.streams.(STREAM_STORE2).filtered');
F465 = zeros(size(allSignals(:,1:N:end-N+1)));
for ii = 1:size(allSignals,1)
    F465(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
end
minLength2 = size(F465,2);

% Create mean signal, standard error of signal, and DC offset of 465 signal
meanSignal2 = mean(F465);
stdSignal2 = std(double(F465))/sqrt(size(F465,1));
dcSignal2 = mean(meanSignal2);

%% Plot Epoch Averaged Response

% Create the time vector for each stream store
ts1 = TRANGE(1) + (1:minLength1) / data.streams.(STREAM_STORE1).fs*N;
ts2 = TRANGE(1) + (1:minLength2) / data.streams.(STREAM_STORE2).fs*N;

% Subtract DC offset to get signals on top of one another
meanSignal1 = meanSignal1 - dcSignal1;
meanSignal2 = meanSignal2 - dcSignal2;

% Plot the 405 and 465 average signals
figure(1);
% subplot(4,1,1)
plot(ts1, meanSignal1, 'color',[0.4660, 0.6740, 0.1880], 'LineWidth', 3); hold on;
plot(ts2, meanSignal2, 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); 


% Plot vertical line at epoch onset, time = 0
line([0 0], [min(F465(:) - dcSignal2), max(F465(:)) - dcSignal2], 'Color', [.7 .7 .7], 'LineStyle','-', 'LineWidth', 3)

% Make a legend
legend('405 nm','465 nm','Shock Onset', 'AutoUpdate', 'off');

% Create the standard error bands for the 405 signal
XX = [ts1, fliplr(ts1)];
YY = [meanSignal1 + stdSignal1, fliplr(meanSignal1 - stdSignal1)];

% Plot filled standard error bands.
h = fill(XX, YY, 'g');
set(h, 'facealpha',.25,'edgecolor','none')

% Repeat for 465
XX = [ts2, fliplr(ts2)];
YY = [meanSignal2 + stdSignal2, fliplr(meanSignal2 - stdSignal2)];
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')

% Finish up the plot
axis tight
xlabel('Time, s','FontSize',12)
ylabel('V', 'FontSize', 12)

title(sprintf('Big Reward Session 3_from Filter, %d Trials (%d Artifacts Removed)', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts))
set(gcf, 'Position',[100, 100, 800, 500])

hold off;

%%
%CALC 

[numTrials,~]=size(BehavData.collectionTime(:));
Tris=[1:numTrials]';


% if EPOC == 'Collect';
%     BL_time = BehavData.collectionTime(:)-BehavData.choiceTime(:);
% end


%%

% Heat Map based on z score of 405 fit subtracted 465

% Fitting 405 channel onto 465 channel to detrend signal bleaching
% Scale and fit data
% Algorithm sourced from Tom Davidson's Github:
% https://github.com/tjd2002/tjd-shared-code/blob/master/matlab/photometry/FP_normalize.m




bls = polyfit(F405(1:end), F465(1:end), 1); %polyfit(F465(1:end), F405(1:end), 1);
Y_fit_all = bls(1) .* F405 + bls(2);
Y_dF_all = F465 - Y_fit_all;

zall = zeros(size(Y_dF_all));
tmp = 0;


% for i = 1:size(Y_dF_all,1)
%     BL_shifted=[BASELINE_PER(1)+BL_time(i) BASELINE_PER(2)+BL_time(i)];
%     ind = ts2(1,:) < BL_shifted(2) & ts2(1,:) > BL_shifted(1);
%     zb = mean(Y_dF_all(i,ind)); % baseline period mean (-10sec to -6sec)
%     zsd = std(Y_dF_all(i,ind)); % baseline period stdev
%     for j = 1:size(Y_dF_all,2) % Z score per bin
%         tmp = tmp + 1;
%         zall(i,tmp)=(Y_dF_all(i,j) - zb)/zsd;
%     end
%     tmp=0;
% end

pp = 1;

% Comment out BL_shifted(pp,:) and the associated ind to not adjust
% baseline to the ITI period

for i = 1:size(Y_dF_all,1)
%     BL_shifted(pp,:)=[BASELINE_PER(1)+time2EPOC(i) BASELINE_PER(2)+time2EPOC(i)]; %BL_shifted(pp,:)=[BASELINE_PER(1)+(-1*time2Collect(i)) BASELINE_PER(2)+time2Collect(i)];
%     ind = ts2(1,:) < BL_shifted(pp,2) & ts2(1,:) > BL_shifted(pp,1);
    ind = ts2(1,:) < BASELINE_PER(2) & ts2(1,:) > BASELINE_PER(1);
    zb = mean(Y_dF_all(i,ind)); % baseline period mean
    zbmedian = median(Y_dF_all(i,length(ts1)));
    zsd = std(Y_dF_all(i,ind)); % baseline period stdev
    pp=pp+1;
    for j = 1:size(Y_dF_all,2) % Z score per bin
        tmp = tmp + 1;
        zall(i,tmp)=(Y_dF_all(i,j) - zb)/zsd;
        dfALL(i,tmp)=(Y_dF_all(i,j) - zbmedian)/zbmedian;
        BL_mean(i) = mean(Y_dF_all(i,:));
    end
    tmp=0;
end
%%
% Standard error of the z-score
zerror = std(zall)/size(zall,1);
ZallMean=mean(zall,1);
dfALLmean = mean(dfALL,1);
% Plot heat map
% subplot(4,1,2)
figure(2);
IM=imagesc(ts2, 1, zall);hold on;

load('tokyo.mat'); %using Scientific Colour-Maps 6.0 (http://www.fabiocrameri.ch/colourmaps.php)
colormap(tokyo); % c1 = colorbar; 
scatter(time2Collect,Tris,'Marker','p','MarkerFaceColor','w')
plot(zeros(numTrials,1),Tris)
%scatter(time2choose,Tris,'Marker','>','MarkerFaceColor','k')
title(sprintf('Z-Score Heat Map, %d Trials (%d Artifacts Removed)', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts));
ylabel('Trials', 'FontSize', 12);
hold off;


%%
figure(3);
% Fill band values for second subplot. Doing here to scale onset bar
% correctly
XX = [ts2, fliplr(ts2)];
YY = [mean(zall)-zerror, fliplr(mean(zall)+zerror)];

% subplot(4,1,3)

plot(ts2, mean(zall), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); hold on;
line([0 0], [min(YY), max(YY)], 'Color', [.7 .7 .7], 'LineWidth', 2)

h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')

% Finish up the plot
axis tight
xlabel('Time, s','FontSize',12)
ylabel('Z-score', 'FontSize', 12)
title(sprintf('465 nm Foot Shock Response, %d Trials (%d Artifacts Removed)', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts))
%c2 = colorbar;

hold off;

%% Plot z-scored traces for each trial
figure(4);
plot(ts1,zall)

%%
% % Quantify changes as area under the curve for cue (-5 sec) and shock (0 sec)
AUC=[]; % cue, shock
AUC(1,1)=trapz(mean(zall(:,ts2(1,:) < -0 & ts2(1,:) > -2)));
AUC(1,2)=trapz(mean(zall(:,ts2(1,:) > 0 & ts2(1,:) < 1)));
subplot(4,1,4);
hBar = bar(AUC, 'FaceColor', [.8 .8 .8]);

% Run a two-sample T-Test
[h,p,ci,stats] = ttest2(mean(zall(:,ts2(1,:) < -3 & ts2(1,:) > -5)),mean(zall(:,ts2(1,:) > 0 & ts2(1,:) < 2)));

% Plot significance bar if p < .05
hold on;
centers = get(hBar, 'XData');
plot(centers(1:2), [1 1]*AUC(1,2)*1.1, '-k', 'LineWidth', 2)
p1 = plot(mean(centers(1:2)), AUC(1,2)*1.2, '*k');
set(gca,'xticklabel',{'Cue','Shock'});
title({'Cue vs Shock Response Changes', 'Area under curve'})
legend(p1, 'Check Statistics','Location','southeast');

set(gcf, 'Position',[100, 100, 500, 1000])
%%

c = distinguishable_colors(numel(BehavData.Trial));

figure(5)
set(gca, 'ColorOrder', c, 'NextPlot', 'replacechildren');

plot(ts1, zall, 'LineWidth',1);

ID_legend = num2str(BehavData.Trial);
legend(ID_legend);