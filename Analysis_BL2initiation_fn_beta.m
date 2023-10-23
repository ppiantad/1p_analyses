function [data, ZallMean,zerror,ts1,numTrials,STREAM_STORE1,STREAM_STORE2,TTLNumber, Descriptives, BehavData, zall] =Analysis_BL2initiation_fn_beta(animalName,blockpath,whichStreams,whichTTL,behavFiles, EPOC);
 


%%

data = TDTbin2mat(blockpath);
%%

%work behavioral data into photo data, generate descriptive statistics
[BehavData,ABETfile, Descriptives]=ABET2TableFn_Chamber_Av3(behavFiles);

%align timestamps of behavioral data to timestamps of photometry
%find time of first TTL.  This value will be passed to Behavioral Data so
%that the timestamps of the ABET data will match the timestamps of Synapse

TTLNumber=fieldnames(data.epocs);
if whichTTL==2
   epoc1=2; 
elseif whichTTL==1
   epoc1=1;
elseif whichTTL==3
   epoc1=3;
elseif whichTTL==4
    epoc1=4;
elseif whichTTL==5
    epoc1=5;
end

timeStart=data.epocs.(TTLNumber{epoc1}).onset(1);
timeShift=timeStart*ones(numel(BehavData.choiceTime(:)),1);
BehavData.choiceTime(:)=BehavData.choiceTime(:)+timeShift;
BehavData.collectionTime(:)=BehavData.collectionTime(:)+timeShift;
BehavData.stTime(:)=BehavData.stTime(:)+timeShift;



%%
% TTLNumber=fieldnames(data.epocs);
% if whichTTL==2
%     epoc1=1; 
% elseif whichTTL==1
%    epoc1=1;
% end





%%
%pick trial types
% If analyzing entire dataset for a given exp, set to 'ALL', 1
BehavData=TrialFilter(BehavData,'ALL',1);


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
data.epocs.start.data = [1:numel(BehavData.stTime)]';



%%

%find channel names
streamNames=fieldnames(data.streams);
if whichStreams==12
    stream1=2; stream2=1;
elseif whichStreams==34
    stream1=3; stream2=4;
end

%%
%use TDT filter to get just the trials you want
% EPOC = 'Choice'; %set which epoc to look at THIS WAS WRONG BEFORE
STREAM_STORE1 = streamNames{stream1}; %put here what your 405 channel is
STREAM_STORE2 = streamNames{stream2}; %put here what your 465 channel is
TRANGE = [-33 43]; % window size [start time relative to epoc onset, window duration](if baselining to START, which is initiation) has to include enough time to account for 30s pre-choice time + baseline period)  -10 20 default -4 12 -31 41 -41 51 -51 61 -33 43 
BASELINE_PER = [-10 -5]; % baseline period within our window -6 -1 default; THE FIRST VALUE CAN'T BE < THE FIRST TRANGE VALUE -10 -5 -2 -0 -3 -1 -3 0
ARTIFACT = Inf; % optionally set an artifact rejection level
data = TDTfilter(data,EPOC,'TIME',TRANGE); %perform TDTfilter 

%%
time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);

if strcmpi(EPOC, 'start')
    time2EPOC = BehavData.stTime(:) - BehavData.choiceTime(:);
elseif strcmpi(EPOC, 'Collect')
    BL_time = BehavData.collectionTime(:)-BehavData.choiceTime(:);
    time2EPOC = BehavData.choiceTime(:)-BehavData.collectionTime(:); %time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);
elseif strcmpi(EPOC, 'Choice')
    BL_time = BehavData.stTime(:)-BehavData.choiceTime(:);
%     time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);
    time2EPOC = BehavData.stTime(:)-BehavData.choiceTime(:);
end


%%

%Check if any trials exist

if isempty(BehavData)
    ZallMean=999;
    zerror=999;
    ts1=999;
    numTrials=999; %set numTrials to 999 because this is an impossible # of trials, can then be filtered out for behavioral analysis in other code
    zall=999;
% if isempty(data.epocs.(EPOC).data)
%     ZallMean=999;
%     zerror=999;
%     ts1=999;
%     numTrials=0;
   
else


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

   % Plot Epoch Averaged Response

    % Create the time vector for each stream store
    ts1 = TRANGE(1) + (1:minLength1) / data.streams.(STREAM_STORE1).fs*N;
    ts2 = TRANGE(1) + (1:minLength2) / data.streams.(STREAM_STORE2).fs*N;

    % Subtract DC offset to get signals on top of one another
    meanSignal1 = meanSignal1 - dcSignal1;
    meanSignal2 = meanSignal2 - dcSignal2;


%ADDED TO BL TO NSP 

[numTrials,~]=size(BehavData.collectionTime(:));
Tris=[1:numTrials]';

% time2Collect=BehavData.collectionTime(:)-BehavData.choiceTime(:);
% time2choose=BehavData.stTime(:)-BehavData.choiceTime(:);

    % Heat Map based on z score of 405 fit subtracted 465

    % Fitting 405 channel onto 465 channel to detrend signal bleaching
    % Scale and fit data
    % Algorithm sourced from Tom Davidson's Github:
    % https://github.com/tjd2002/tjd-shared-code/blob/master/matlab/photometry/FP_normalize.m

    bls = polyfit( F405(1:end),F465(1:end), 1);
    Y_fit_all = bls(1) .* F405 + bls(2);
    Y_dF_all = F465 - Y_fit_all;
    
    zall = zeros(size(Y_dF_all));
    tmp = 0;
    

for i = 1:size(Y_dF_all,1)
%     BL_shifted=[BASELINE_PER(1)+time2EPOC(i) BASELINE_PER(2)+time2EPOC(i)]; %BL_shifted(pp,:)=[BASELINE_PER(1)+(-1*time2Collect(i)) BASELINE_PER(2)+time2Collect(i)];
    ind = ts2(1,:) < BASELINE_PER(2) & ts2(1,:) > BASELINE_PER(1); %use if you want to test without shifting baseline according to particular EPOC
%     ind = ts2(1,:) < BL_shifted(2) & ts2(1,:) > BL_shifted(1);
    zb = mean(Y_dF_all(i,ind)); % baseline period mean
%     zbmedian = median(Y_dF_all(i,length(ts1)));
    zsd = std(Y_dF_all(i,ind)); % baseline period stdev
    for j = 1:size(Y_dF_all,2) % Z score per bin
        tmp = tmp + 1;
        zall(i,tmp)=(Y_dF_all(i,j) - zb)/zsd;
%         dfALL(i,tmp)=(Y_dF_all(i,j) - zbmedian)/zbmedian;
%         BL_mean(i) = mean(Y_dF_all(i,:));

    end
    tmp=0;
end

    % Standard error of the z-score
    zerror = std(zall)/size(zall,1);
    
    
    ZallMean=mean(zall,1);
    [numTrials, ~] = size(F465);
end   
end
