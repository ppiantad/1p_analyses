% need to run first section of patDataScript.m first, with TrialFilter set
% to ALL
% plus need to run SLEAP_attempts_inscopix.m
animalIDs = (fieldnames(final));


% String to compare
targetAnimal = 'BLA_Insc_24';

% Perform element-wise comparison
animal_index_to_plot = find(strcmp(animalIDs, targetAnimal));

session_to_analyze = 'RDT_D1';

ca_data_type = "C_raw"; % C % C_raw %S
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes

epoc_to_align = 'choiceTime';

neuron_num_to_model = 27;

close all
%%
ca = final.(targetAnimal).(session_to_analyze).CNMFe_data.(ca_data_type);
zscored_SLEAP_data_vel_filtered_session = final_SLEAP.(targetAnimal).(session_to_analyze).zscored_SLEAP_data_velocity;
velocity_time = final_SLEAP.(targetAnimal).(session_to_analyze).SLEAP_data.idx_time;

uv = final.(targetAnimal).(session_to_analyze).uv;
uv.evtWin = [-8 8]; %what time do you want to look at around each event [-2 8] [-10 5]
uv.BLper = [-10 -5];
BehavData =uv.BehavData;
BehavData = BehavData(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3, :);

frames3 = final.(targetAnimal).(session_to_analyze).time  ;
ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);


%preloads basis set
load('basis_81x25.mat')

sampling_rate = uv.dt*100;

%how much time should you shift back (in seconds)
time_back_orig=10; % 2 original
time_forward_orig=10; % 6 original

time_window = (time_back_orig+time_forward_orig)*sampling_rate; 


% Specify the number of rows you want (200 in this case)
desired_rows = 201;

% Replicate the basis set to create a 200x25 basis set
extended_basis_set = repmat(basis_set, ceil(desired_rows / size(basis_set, 1)), 1);

% Trim the excess rows
extended_basis_set = extended_basis_set(1:desired_rows, :);

basis_set = extended_basis_set;

%Z-score gcamp signal PTP: isn't this the processed GCAMP? wtf
% I am currently using this because this is what the code on github from
% the Witten lab (Parker 2022) used! Change to "raw" GCAMP might be more
% legit. 
% gcamp_y=(neuron.C(neuron_num_to_model,:)-mean(neuron.C(neuron_num_to_model,:)))./std(neuron.C(neuron_num_to_model,:));
% gcamp_y=(neuron.C_raw(neuron_num_to_model,:)-mean(neuron.C_raw(neuron_num_to_model,:)))./std(neuron.C_raw(neuron_num_to_model,:));
gcamp_y=(ca(neuron_num_to_model,:)-mean(ca(neuron_num_to_model,:)))./std(ca(neuron_num_to_model,:));
% gcamp_y_filt = sgolayfilt(gcamp_y, 1,75);
% gcamp_y = gcamp_y_filt; 


% trim gcamp_y to length of motion data, or vice versa, depending on when
% behavioral video was stopped relative to Ca recording

if length(gcamp_y) > length(zscored_SLEAP_data_vel_filtered_session)
    gcamp_y = gcamp_y(1,1:length(zscored_SLEAP_data_vel_filtered_session));
elseif length(zscored_SLEAP_data_vel_filtered_session) > length(gcamp_y)
    zscored_SLEAP_data_vel_filtered_session = zscored_SLEAP_data_vel_filtered_session(1,1:length(gcamp_y));
end





%how many shuffled iterations to make distribution of F-stats
num_shuff=500;

%load gcamp and behavioral data of example recording
% load('example_neuron.mat')




%convert seconds to hertz for relevant behavioral timestamps
Large_yes_shk=BehavData.choiceTime(BehavData.shock == 1 & BehavData.bigSmall == 1.2);
Large_no_shk=BehavData.choiceTime(BehavData.shock ~= 1 & BehavData.bigSmall == 1.2);
Small_no_shk = BehavData.choiceTime(BehavData.shock ~= 1 & BehavData.bigSmall == 0.3);


Times_TrialStart = (BehavData.stTime*sampling_rate)'; 
Times_Large_no_shk=(Large_no_shk*sampling_rate)';
Times_Large_yes_shk=(Large_yes_shk*sampling_rate)';
Times_Small_no_shk_redo=(Small_no_shk*sampling_rate);
Times_RewardCollection=(BehavData.collectionTime*sampling_rate)'; 
velocity = zscored_SLEAP_data_vel_filtered_session;



%define time shift for each tested event, 0 for -2-6 seconds, 1 for 0-8 seconds
con_shift=[0 0 0 0 0 0 0];

figure;
plot(velocity_time, velocity); hold on;
plot(velocity_time, gcamp_y)

% 
% LeverPresent=output.LeverPresentation.*g_output.samp_rate;
% LeverTimesI=output.LeverTimes(output.IpsPress==1).*g_output.samp_rate;
% LeverTimesC=output.LeverTimes(output.IpsPress==-1).*g_output.samp_rate;
% CSRew=output.RewardPresentation(output.RewardEnter~=0).*g_output.samp_rate;
% CSNoRew=output.RewardPresentation(output.RewardEnter==0).*g_output.samp_rate;
% RewardEnter=output.RewardEnter(output.RewardEnter~=0).*g_output.samp_rate;

%Define behaviors to be tested in model
cons={'Times_TrialStart','Times_Large_no_shk','Times_Large_yes_shk','Times_Small_no_shk_redo'...
    'Times_RewardCollection', 'velocity'};

%% Regression data prep

%Initialize x-matrices
con_iden=[];
x_basic=[];    %No interaction terms, simply event times
event_times_mat=[];
num_bins=numel(gcamp_y);

for con=1:numel(cons)

    if con_shift(con)==1
        time_back=0;
        time_forward=time_back_orig+time_forward_orig;
    else
        time_back=time_back_orig;
        time_forward=time_forward_orig;
    end

    %gets matrix of event times (in hertz)
    con_times=eval(cons{con});

    %Gets rid of abandonded trials
    con_times(con_times==0)=[];

    %Creates vector with binary indication of events
    con_binned=zeros(1,num_bins);
    if ~strcmp(cons{con}, 'velocity')
        con_binned(int32(con_times))=1;
        % PTP: this appears to shift the index for each event back to the
        % appropriate time_back spot
        % con_binned=circshift(con_binned,[0,-time_back*g_output.samp_rate]);
        con_binned=circshift(con_binned,[0,-time_back*sampling_rate]);
        event_times_mat=vertcat(event_times_mat,con_binned);
    elseif strcmp(cons{con}, 'velocity')
        con_binned = velocity;
    end




    %convolves time series of events with basis sets
    % PTP these are vertically concatenated (each convolved set is 25 rows
    % wide)
    if strcmp(cons{con}, 'velocity')
        x_basic=horzcat(x_basic,velocity');
    elseif ~strcmp(cons{con}, 'velocity')
        for num_sets=1:numel(basis_set(1,:))
            temp_conv_vec=conv(con_binned,basis_set(:,num_sets));
            x_basic=horzcat(x_basic,temp_conv_vec(1:numel(con_binned))');


        end
    end

    if strcmp(cons{con}, 'velocity')
        con_iden = [con_iden con];

    elseif ~strcmp(cons{con}, 'velocity')
        con_iden=[con_iden ones(1,size(basis_set,2))*con];
    end
end

%mean center predictor matrix
x_all=mean_center(x_basic);

%makes predictor cell, identifier of groups for nested regressions
num_cons=unique(con_iden);
for i=1:numel(num_cons)
    preds_cells{i}=find(con_iden==i);
end


%% First, unshuffled F-stat
X = x_all;
y = gcamp_y';
preds_to_test_cell = preds_cells;

X = [ones(size(X,1),1) X];
n = size(X,1);  % number of observations
k = size(X,2);             % number of variables (including constant)

% Perform Lasso regression
[B, FitInfo] = lasso(X, y, 'Alpha', 1);

% Find the lambda that minimizes cross-validated MSE
[~, idxLambdaMinMSE] = min(FitInfo.MSE);

% Extract the coefficients for the lambda that minimizes MSE
b_lasso = B(:, idxLambdaMinMSE);

% Calculate residuals
u = y - X * b_lasso;

% Calculate F-statistics
clear Fp_vec
for l=1:length(preds_to_test_cell)
    preds_to_test = preds_to_test_cell{l} + 1; % because of the constant

    R = zeros(length(preds_to_test), size(X, 2));
    for l2 = 1:length(preds_to_test)
        R(l2, preds_to_test(l2)) = 1;
    end

    r = zeros(length(preds_to_test), 1);  % Testing restriction: R * b = r

    num_restrictions = size(R, 1);
    F = (R * b_lasso - r)' * inv(R * (X' * X) * R') * (R * b_lasso - r) / num_restrictions;  % F-stat (see Hiyashi for reference)
    F_vec(l) = F;
    Fp_vec(l) = 1 - fcdf(F, num_restrictions, size(X, 1) - size(X, 2));  % F p-val
end

Fp_nonshuff = Fp_vec;
F_nonshuff = F_vec;
b_final= b_lasso';

%% Next, shuffled F-stats
for i=1:num_shuff
    X = x_all;
    y=circshift(gcamp_y',[randi(numel(gcamp_y'),1)]);
    preds_to_test_cell = preds_cells;

    X = [ones(size(X,1),1) X];
    n = size(X,1);  % number of observations
    k = size(X,2);             % number of variables (including constant)

    b = X \ y;                 % estimate b with least squares
    u = y - X * b;             % calculates residuals
    s2 = u' * u / (n - k);     % estimate variance of error term (assuming homoskedasticity, independent observations)
    BCOV = inv(X'*X) * s2;     % get covariance matrix of b assuming homoskedasticity of error term etc...
    bse = diag(BCOV).^.5;      % standard errors

    clear Fp_vec
    for l=1:length(preds_to_test_cell)
        preds_to_test_cell{l} = preds_to_test_cell{l}+1; %because of the constant

        R = zeros(length(preds_to_test_cell{l}),k);
        for l2 = 1:length(preds_to_test_cell{l})
            R(l2, preds_to_test_cell{l}(l2))=1;
        end

        r = zeros(length(preds_to_test_cell{l}),1);          % Testing restriction: R * b = r

        num_restrictions = size(R, 1);
        F = (R*b - r)'*inv(R * BCOV * R')*(R*b - r) / num_restrictions;   % F-stat (see Hiyashi for reference)
        F_vec(l) = F;
        Fp_vec(l) = 1 - fcdf(F, num_restrictions, n - k);  % F p-val
    end

    Fp_shuff(:,i)=Fp_vec;
    F_shuff(:,i)=F_vec;

end

%% Compute p-values from shuffled distribution of F-stats
for g=1:numel(preds_cells)
    temp_sort=squeeze(sort(F_shuff(g,:)));
    temp_pval=find(F_nonshuff(g)<temp_sort); if isempty(temp_pval)==1 ;temp_pval=numel(temp_sort) ;end
    pvals(g)=1-temp_pval(1)/numel(temp_sort);
end

%% Output structure
encoding_output.coefs = b_final; %coefficients for each event 
encoding_output.pvals = pvals; %p-values for each event
encoding_output.events = cons; %behaviors tested
encoding_output.event_identiy = preds_cells; %cell of identifiers for the coefficients of each event


%%
figure; plot(encoding_output.coefs(1, 2:26));
figure; plot(encoding_output.coefs(1, 27:51));
figure; plot(encoding_output.coefs(1, 52:76));
figure; plot(encoding_output.coefs(1, 77:101));
figure; plot(encoding_output.coefs(1, 102:126));



first = encoding_output.coefs(1, 2:26);
second = encoding_output.coefs(1, 27:51);
third = encoding_output.coefs(1, 52:76);
fourth = encoding_output.coefs(1, 77:101);
fifth = encoding_output.coefs(1, 102:126);
all_coefs = [first; second; third; fourth; fifth];
sum_coeffs = sum(all_coefs);

figure; plot(sum_coeffs);

%%
unitTrace = ca(neuron_num_to_model,:); %get trace
eTS = BehavData.(epoc_to_align); %get time stamps
for t = 1:size(eTS,1)
    %% set each trial's temporal boundaries
    timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
    BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
    if min(timeWin) > min(frames3) & max(timeWin) < max(frames3)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
        %% get unit event counts in trials
        %% get unit ca traces in trials
        idx = frames3 > min(timeWin) & frames3 < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
        sum(idx)
        bl_idx = frames3 > min(BL_win) & frames3 < max(BL_win);
        %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
        caTraceTrials(t,1:sum(idx)) = unitTrace(idx);
        % zb(t,:) = mean(unitTrace(bl_idx)); %baseline mean
        % zsd(t,:) = std(unitTrace(bl_idx)); %baseline std
        zb(t,:) = mean(caTraceTrials(t,:)); %baseline mean
        zsd(t,:) = std(caTraceTrials(t,:)); %baseline std


        tmp = 0;
        for j = 1:size(caTraceTrials,2)
            tmp = tmp+1;
            zall(t,tmp) = (caTraceTrials(t,j) - zb(t))/zsd(t);
        end
        clear j;



    end
end
zall = zall(:, 1:size(ts1, 2));
ca_mean = mean(zall); 

figure; plot(ca_mean(1, 1:160)'); 
figure;
imagesc(ts1, 1, zall);
title('Z-scored Trial Data (normalized)')