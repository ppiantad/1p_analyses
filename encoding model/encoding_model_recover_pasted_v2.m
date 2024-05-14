% load('BLA_Insc_24_SLEAP_data.mat')
% load('BLA_Insc_24_Ca_data.mat')

%%

animalIDs = (fieldnames(final));

% String to compare
targetAnimal = 'BLA_Insc_24';

session_to_analyze = 'RDT_D1';

ca_data_type = "C_raw"; % C % C_raw %S
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes
% CNMFe_data.spike_prob: CASCADE inferred spikes - multiply x sampling rate
% (10) for spike rate

epoc_to_align = 'choiceTime';

neuron_num_to_model = 2; %50

ca = final.(targetAnimal).(session_to_analyze).CNMFe_data.(ca_data_type);
zscored_SLEAP_data_vel_filtered_session = final_SLEAP.(targetAnimal).(session_to_analyze).zscored_SLEAP_data_velocity;
velocity_time = final_SLEAP.(targetAnimal).(session_to_analyze).SLEAP_data.idx_time;

uv = final.(targetAnimal).(session_to_analyze).uv;
uv.evtWin = [-8 8]; %what time do you want to look at around each event [-2 8] [-10 5]
uv.BLper = [-10 -5];
BehavData = final.(targetAnimal).(session_to_analyze).uv.BehavData;
% BehavData = BehavData(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3, :);

% Check if 'type_binary' exists in the column headers
isTypeBinary = ismember('type_binary', BehavData.Properties.VariableNames);

if ~isTypeBinary
    % Add a new column with header 'type_binary' containing some default values
    defaultValues = NaN(height(BehavData), 1); % Change NaN to whatever default value you need
    BehavData = addvars(BehavData, defaultValues, 'NewVariableNames', 'type_binary');
    disp('The column header "type_binary" was added to the table.');
else
    disp('The column header "type_binary" already exists in the table.');
end

frames3 = final.(targetAnimal).(session_to_analyze).time;
ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);

% %preloads basis set
% load('basis_81x25.mat')
%preloads basis set
load('updated_basis_set.mat')

sampling_rate = uv.dt*100;

%how much time should you shift back (in seconds)
time_back_orig=10; % 2 original
time_forward_orig=10; % 6 original

time_window = (time_back_orig+time_forward_orig)*sampling_rate; 

% Z-score gcamp signal PTP: isn't this the processed GCAMP? wtf
% I am currently using this because this is what the code on github from
% the Witten lab (Parker 2022) used! Change to "raw" GCAMP might be more
% legit. 
% gcamp_y=(neuron.C(neuron_num_to_model,:)-mean(neuron.C(neuron_num_to_model,:)))./std(neuron.C(neuron_num_to_model,:));
% gcamp_y=(neuron.C_raw(neuron_num_to_model,:)-mean(neuron.C_raw(neuron_num_to_model,:)))./std(neuron.C_raw(neuron_num_to_model,:));
gcamp_y=(ca(neuron_num_to_model,:)-nanmean(ca(neuron_num_to_model,:)))./nanstd(ca(neuron_num_to_model,:));
% gcamp_y_zscored = zscore(ca(neuron_num_to_model,:));
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
num_shuff=500; %500

%convert seconds to frames for relevant behavioral timestamps
Large_yes_shk=BehavData.choiceTime(BehavData.shock == 1 & BehavData.bigSmall == 1.2);
Large_no_shk=BehavData.choiceTime(BehavData.shock ~= 1 & BehavData.bigSmall == 1.2);
Small_no_shk = BehavData.choiceTime(BehavData.shock ~= 1 & BehavData.bigSmall == 0.3);
valid_rew_collect = BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3);
valid_trial_starts = BehavData.stTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3 | BehavData.omissionALL == 1);

blank_touches = BehavData.choiceTime(BehavData.Blank_Touch == 1);
large_aborts = BehavData.choiceTime(BehavData.type_binary == 1);
small_aborts = BehavData.choiceTime(BehavData.type_binary == 2);
lose_shifts = BehavData.choiceTime(BehavData.lose_shift == 1);
omissions = BehavData.choiceTime(BehavData.omissionALL == 1);

Times_TrialStart = (valid_trial_starts*sampling_rate)'; 
Times_Large_no_shk=(Large_no_shk*sampling_rate)';
Times_Large_yes_shk=(Large_yes_shk*sampling_rate)';
Times_Small_no_shk_redo=(Small_no_shk*sampling_rate)';
Times_RewardCollection=(valid_rew_collect*sampling_rate)'; 
Times_blank_touch = (blank_touches*sampling_rate)';
Times_large_aborts = (large_aborts*sampling_rate)';
Times_small_aborts = (small_aborts*sampling_rate)';
% Times_lose_shifts = (lose_shifts*sampling_rate)';
Times_omissions = (omissions*sampling_rate)';
velocity = zscored_SLEAP_data_vel_filtered_session;

%define time shift for each tested event, 0 for -2-6 seconds, 1 for 0-8 seconds
con_shift=[0 0 0 0 0 0 0 0 0 0];

figure;
plot(velocity_time, velocity); hold on;
plot(velocity_time, gcamp_y)

% LeverPresent=output.LeverPresentation.*g_output.samp_rate;
% LeverTimesI=output.LeverTimes(output.IpsPress==1).*g_output.samp_rate;
% LeverTimesC=output.LeverTimes(output.IpsPress==-1).*g_output.samp_rate;
% CSRew=output.RewardPresentation(output.RewardEnter~=0).*g_output.samp_rate;
% CSNoRew=output.RewardPresentation(output.RewardEnter==0).*g_output.samp_rate;
% RewardEnter=output.RewardEnter(output.RewardEnter~=0).*g_output.samp_rate;

%Define behaviors to be tested in model
cons={'Times_TrialStart','Times_Large_no_shk','Times_Large_yes_shk','Times_Small_no_shk_redo'...
    'Times_RewardCollection', 'Times_blank_touch', 'Times_large_aborts', 'Times_small_aborts', 'Times_omissions', 'velocity'};

eTS = BehavData.(epoc_to_align); %get time stamps

% if you want to check a specific event, uncomment below and re-run
% eTS = valid_rew_collect; 

unitTrace = ca(neuron_num_to_model,:); %get trace

%uncomment below if you want to make a heatmap of a kernel (first_kernel is
%just a random one, can add more below)
% unitTrace = first_kernel'; %get trace
for t = 1:size(eTS,1)
    % set each trial's temporal boundaries
    timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
    BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
    if min(timeWin) > min(frames3) & max(timeWin) < max(frames3)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
        % get unit event counts in trials
        % get unit ca traces in trials
        idx = frames3 > min(timeWin) & frames3 < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
        sum(idx);
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
            % zall(t, tmp) = caTraceTrials(t,j);
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
clear zall ca_mean caTraceTrials

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

Fp_nonshuff(:)=Fp_vec;
F_nonshuff(:)=F_vec;
b_final(:)=b;

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

%the below function is what Ruairi prefers to correct for multiple
%comparisons
adjusted_p_values = mafdr(pvals, 'BHFDR', 'true');
encoding_table = array2table([pvals; adjusted_p_values], 'VariableNames', cons, 'RowNames', {'uncorrected pvals', 'corrected pvals'});

%%
% smooth coefs
encoding_output.coefs = sgolayfilt(encoding_output.coefs(1, :), 9, 21);
begin_val = 2;
for preds = 1:size(preds_cells, 2)
    length_current_pred = length(preds_cells{preds});
  
    figure;
    plot(encoding_output.coefs(1, begin_val:(length_current_pred*preds)+1));
    coefs_by_predictors{preds} = encoding_output.coefs(1, begin_val:(length_current_pred*preds)+1);
    % title(sprintf('Coefs %d to Coefs %d\n', begin_val, (length_current_pred*preds)+1));
    begin_val = begin_val + length(preds_cells{preds});
    title(strrep(cons(preds), '_', ' '));
end


% empty matrix to store the concatenated coefficients
concatenated_coefs = [];


for preds = 1:numel(coefs_by_predictors)
    % Check if the cell is empty
    if ~isempty(coefs_by_predictors{preds})
        % Concatenate the contents horizontally
        concatenated_coefs = [concatenated_coefs; coefs_by_predictors{preds}];
    end
end

sum_coeffs = sum(concatenated_coefs);

figure; plot(sum_coeffs); title('overall coefs')

%%
scaled_coefs = b_final(:,2:end).*x_basic;

all_kernel = sum(scaled_coefs, 2);
figure; plot(gcamp_y);
hold on; plot(all_kernel);
full_kernel_corr_gcamp = corrcoef(all_kernel', gcamp_y)


scaled_coefs_means = mean(scaled_coefs);


first_kernel = mean(scaled_coefs(:, 101:125), 2);

figure; plot(first_kernel)
xline(Times_RewardCollection)


begin_val = 1;
for preds = 1:size(preds_cells, 2)
    length_current_pred = length(preds_cells{preds});
  
    figure;
    plot(scaled_coefs_means(1, begin_val:(length_current_pred*preds)));
    scaled_coefs_by_predictors{preds} = scaled_coefs_means(1, begin_val:(length_current_pred*preds));
    % title(sprintf('Coefs %d to Coefs %d\n', begin_val, (length_current_pred*preds)+1));
    begin_val = begin_val + length(preds_cells{preds});
    title(strrep(cons(preds), '_', ' '));
end


% empty matrix to store the concatenated coefficients
concatenated_scaled_coefs = [];


for preds = 1:numel(coefs_by_predictors)
    % Check if the cell is empty
    if ~isempty(coefs_by_predictors{preds})
        % Concatenate the contents horizontally
        concatenated_scaled_coefs = [concatenated_scaled_coefs; scaled_coefs_by_predictors{preds}];
    end
end

sum_coeffs = sum(concatenated_scaled_coefs);

figure; plot(sum_coeffs); title('overall kernel')

