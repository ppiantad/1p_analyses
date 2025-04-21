
% 
% num_samples = size(final.b78764.D1_Afternoon.CNMFe_data.C_raw, 2);
% 
% time_step = final.b78764.D1_Afternoon.uv.sampling_rate/100;
% 
% time_series = (1:num_samples) * time_step;
% 
% figure; plot(time_series, final.b78764.D1_Afternoon.CNMFe_data.C_raw(1, :))
%%
iter = 0

%%
% load('BLA-NAcShell_Risk_2024_01_04.mat')

%% Edit these uservariables with what you want to look at
uv.evtWin = [-8 8]; %what time do you want to look at around each event [-2 8] [-10 5] [-10 10]
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
% uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at

ca_data_type = "C"; % C % C_raw %S
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes
% CNMFe_data.spike_prob: CASCADE inferred spikes - multiply x sampling rate
% (10) for spike rate


session_to_analyze = 'D3';

% Parameters
session_duration = 12 * 60; % seconds
sampling_rate = 10; % Hz
total_time_points = session_duration * sampling_rate;



recorded_fps = 30;



% Parameters
stimulus_duration = 2 * 60; % 2 minutes in seconds
num_repeats = 3;
total_stimuli = 2;

% Initialize variables
stimulus_times = cell(total_stimuli, 1);
current_time = 0;

session_length_in_min = 12;

% Loop through each stimulus alternately and calculate start and end times
for j = 1:num_repeats
    for i = 1:total_stimuli
        start_time = current_time;
        end_time = start_time + stimulus_duration;
        
        if isempty(stimulus_times{i})
            stimulus_times{i} = [start_time, end_time];
        else
            stimulus_times{i} = [stimulus_times{i}; start_time, end_time];
        end
        
        current_time = end_time;
    end
end


% Multiply every value in stimulus_times by the FPS
for i = 1:total_stimuli
    stimulus_frames{i} = stimulus_times{i} * sampling_rate;
end




shock_start_time = 4 * 60; % First shock in seconds
shock_interval = 60; % Interval between shocks in seconds
shock_duration = 2; % Duration of each shock in seconds
num_shocks = 6;

% Initialize the footshock array
footshock = zeros(1, total_time_points);

% Calculate the indices for each shock
for i = 0:(num_shocks-1)
    shock_start_idx = (shock_start_time + i * shock_interval) * sampling_rate + 1;
    shock_end_idx = shock_start_idx + shock_duration * sampling_rate - 1;
    footshock(shock_start_idx:shock_end_idx) = 1;
end

% Parameters
shock_start_time = 4 * 60; % First shock in seconds
shock_interval = 60; % Interval between shocks in seconds
shock_duration = 2; % Duration of each shock in seconds
num_shocks = 6;

% Initialize variables
shk_on = zeros(1, num_shocks);
shk_off = zeros(1, num_shocks);

% Calculate on and off times for each shock
for i = 0:(num_shocks-1)
    shk_on(i+1) = shock_start_time + i * shock_interval; % Shock start time in seconds
    shk_off(i+1) = shk_on(i+1) + shock_duration; % Shock end time in seconds
end

yoke_data = 0; % 1, set to 1 if you want to be prompted to yoke the number of trials analyzed, set to 0 otherwise

epoc_to_align = 'choiceTime';
ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);

neuron_num = 0;
use_normalized_time = 0;

clear neuron_mean neuron_sem neuron_num zall_mean zall_array zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized sem_all zall_mean_all 

if strcmp('RM_D1', session_to_analyze)| strcmp('RDT_D1', session_to_analyze) | strcmp('Pre_RDT_RM', session_to_analyze)
    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_29', 'BLA_Insc_38', 'BLA_Insc_39'};

    for i = 1:length(fieldsToRemove)
        if isfield(final, fieldsToRemove{i})
            final = rmfield(final, fieldsToRemove{i});
        end
    end
elseif strcmp('RDT_D2', session_to_analyze)

    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_39'};

    for i = 1:length(fieldsToRemove)
        if isfield(final, fieldsToRemove{i})
            final = rmfield(final, fieldsToRemove{i});
        end
    end
end
%% FILTER TO GET UN-SHUFFLED DATA
iter = iter+1;
neuron_num = 0;
animalIDs = (fieldnames(final));


% Initialize a structure to store calcium event counts
calcium_events = struct();

% Process each animal
animalIDs = fieldnames(final);
for ii = 1:size(animalIDs,1)
    currentanimal = char(animalIDs(ii));
    
    % Check if animal has the session we want to analyze
    if isfield(final.(currentanimal), session_to_analyze)
        disp(['Processing animal: ' currentanimal]);
        
        % Get calcium data
        ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);
        if strcmp(ca_data_type, 'S')
            ca = full(ca);
        end
        
        % Get number of neurons for this animal
        num_neurons = size(ca, 1);
        
        % Create a structure for this animal if it doesn't exist
        if ~isfield(calcium_events, currentanimal)
            calcium_events.(currentanimal) = struct();
        end
        
        % Initialize arrays to store event counts for each stimulus type
        calcium_events.(currentanimal).stimulus_type_1 = zeros(num_neurons, size(stimulus_frames{1}, 1));
        calcium_events.(currentanimal).stimulus_type_2 = zeros(num_neurons, size(stimulus_frames{2}, 1));
        calcium_events.(currentanimal).p_values = zeros(num_neurons, 1);
        calcium_events.(currentanimal).responsive_neurons = zeros(num_neurons, 1); % 0=not responsive, 1=stim1, 2=stim2
        
        % Process each neuron
        for neuron_idx = 1:num_neurons
            neuron_trace = ca(neuron_idx, :);
            
            % Find peaks in the trace (calcium events)
            [peaks, peak_locs] = findpeaks(neuron_trace, 'MinPeakHeight', mean(neuron_trace) + 1*std(neuron_trace), ...
                                          'MinPeakDistance', 3);
            
            % Count events in each stimulus period
            % First stimulus type (safe)
            for stim_idx = 1:size(stimulus_frames{1}, 1)
                frame_start = round(stimulus_frames{1}(stim_idx, 1));
                frame_end = round(stimulus_frames{1}(stim_idx, 2));
                
                % Count peaks in this time window
                events_in_window = sum(peak_locs >= frame_start & peak_locs <= frame_end);
                
                % Store the count
                calcium_events.(currentanimal).stimulus_type_1(neuron_idx, stim_idx) = events_in_window;
            end
            
            % Second stimulus type
            for stim_idx = 1:size(stimulus_frames{2}, 1)
                frame_start = round(stimulus_frames{2}(stim_idx, 1));
                frame_end = round(stimulus_frames{2}(stim_idx, 2));
                
                % Count peaks in this time window
                events_in_window = sum(peak_locs >= frame_start & peak_locs <= frame_end);
                
                % Store the count
                calcium_events.(currentanimal).stimulus_type_2(neuron_idx, stim_idx) = events_in_window;
            end
            
            % Normalize event counts by the duration of stimulus periods
            % This is important if stimulus_type_1 and stimulus_type_2 have different durations
            stim1_duration = (stimulus_frames{1}(:,2) - stimulus_frames{1}(:,1)) / sampling_rate;
            stim2_duration = (stimulus_frames{2}(:,2) - stimulus_frames{2}(:,1)) / sampling_rate;
            
            norm_stim1_events = calcium_events.(currentanimal).stimulus_type_1(neuron_idx, :) ./ stim1_duration';
            norm_stim2_events = calcium_events.(currentanimal).stimulus_type_2(neuron_idx, :) ./ stim2_duration';
            
            % Perform statistical test to compare responses between stimulus types
            % Using paired t-test (if same number of trials) or unpaired t-test (if different number)
            if length(norm_stim1_events) == length(norm_stim2_events)
                [h, p] = ttest(norm_stim1_events, norm_stim2_events);
            else
                [h, p] = ttest2(norm_stim1_events, norm_stim2_events);
            end
            
            % Store p-value
            calcium_events.(currentanimal).p_values(neuron_idx) = p;
            
            % Determine if neuron is responsive (using alpha = 0.05)
            if p < 0.05
                if mean(norm_stim1_events) > mean(norm_stim2_events)
                    calcium_events.(currentanimal).responsive_neurons(neuron_idx) = 1; % Prefers stimulus 1
                else
                    calcium_events.(currentanimal).responsive_neurons(neuron_idx) = 2; % Prefers stimulus 2
                end
            end
            
            % Store normalized event rates
            calcium_events.(currentanimal).norm_event_rate_stim1(neuron_idx, :) = norm_stim1_events;
            calcium_events.(currentanimal).norm_event_rate_stim2(neuron_idx, :) = norm_stim2_events;
        end
        
        % Calculate summary statistics for this animal
        calcium_events.(currentanimal).mean_events_stim1 = mean(calcium_events.(currentanimal).norm_event_rate_stim1, 2);
        calcium_events.(currentanimal).mean_events_stim2 = mean(calcium_events.(currentanimal).norm_event_rate_stim2, 2);
        
        % Calculate percentage of neurons responsive to each stimulus type
        total_neurons = num_neurons;
        stim1_responsive = sum(calcium_events.(currentanimal).responsive_neurons == 1);
        stim2_responsive = sum(calcium_events.(currentanimal).responsive_neurons == 2);
        
        calcium_events.(currentanimal).percent_stim1_responsive = (stim1_responsive / total_neurons) * 100;
        calcium_events.(currentanimal).percent_stim2_responsive = (stim2_responsive / total_neurons) * 100;
        calcium_events.(currentanimal).percent_responsive = ((stim1_responsive + stim2_responsive) / total_neurons) * 100;
        
        disp(['Animal ' currentanimal ': ' num2str(calcium_events.(currentanimal).percent_responsive) '% of neurons show differential responses']);
        disp([' - ' num2str(calcium_events.(currentanimal).percent_stim1_responsive) '% prefer Safe stimulus']);
        disp([' - ' num2str(calcium_events.(currentanimal).percent_stim2_responsive) '% prefer Different stimulus']);
    end
end

% Aggregate results across all animals
all_p_values = [];
all_responsive = [];
all_mean_stim1 = [];
all_mean_stim2 = [];

for ii = 1:length(animalIDs)
    currentanimal = char(animalIDs(ii));
    if isfield(calcium_events, currentanimal)
        all_p_values = [all_p_values; calcium_events.(currentanimal).p_values];
        all_responsive = [all_responsive; calcium_events.(currentanimal).responsive_neurons];
        all_mean_stim1 = [all_mean_stim1; calcium_events.(currentanimal).mean_events_stim1];
        all_mean_stim2 = [all_mean_stim2; calcium_events.(currentanimal).mean_events_stim2];
    end
end

% Apply multiple comparisons correction (False Discovery Rate)
[fdr_p, fdr_mask] = mafdr(all_p_values, 'BHFDR', true);
all_responsive_fdr = zeros(size(all_responsive));
for i = 1:length(fdr_mask)
    if fdr_mask(i)
        if all_mean_stim1(i) > all_mean_stim2(i)
            all_responsive_fdr(i) = 1;
        else
            all_responsive_fdr(i) = 2;
        end
    end
end

% Group summary statistics
group_summary = struct();
group_summary.total_neurons = length(all_p_values);
group_summary.stim1_responsive = sum(all_responsive == 1);
group_summary.stim2_responsive = sum(all_responsive == 2);
group_summary.stim1_responsive_fdr = sum(all_responsive_fdr == 1);
group_summary.stim2_responsive_fdr = sum(all_responsive_fdr == 2);

group_summary.percent_stim1_responsive = (group_summary.stim1_responsive / group_summary.total_neurons) * 100;
group_summary.percent_stim2_responsive = (group_summary.stim2_responsive / group_summary.total_neurons) * 100;
group_summary.percent_responsive = ((group_summary.stim1_responsive + group_summary.stim2_responsive) / group_summary.total_neurons) * 100;

group_summary.percent_stim1_responsive_fdr = (group_summary.stim1_responsive_fdr / group_summary.total_neurons) * 100;
group_summary.percent_stim2_responsive_fdr = (group_summary.stim2_responsive_fdr / group_summary.total_neurons) * 100;
group_summary.percent_responsive_fdr = ((group_summary.stim1_responsive_fdr + group_summary.stim2_responsive_fdr) / group_summary.total_neurons) * 100;

% Display group summary
disp('===== GROUP SUMMARY =====');
disp(['Total neurons analyzed: ' num2str(group_summary.total_neurons)]);
disp(['Uncorrected: ' num2str(group_summary.percent_responsive) '% of neurons show differential responses']);
disp([' - ' num2str(group_summary.percent_stim1_responsive) '% prefer Safe stimulus']);
disp([' - ' num2str(group_summary.percent_stim2_responsive) '% prefer Different stimulus']);
disp(['FDR-corrected: ' num2str(group_summary.percent_responsive_fdr) '% of neurons show differential responses']);
disp([' - ' num2str(group_summary.percent_stim1_responsive_fdr) '% prefer Safe stimulus']);
disp([' - ' num2str(group_summary.percent_stim2_responsive_fdr) '% prefer Different stimulus']);

% Create visualization
figure;

% Plot 1: Proportion of responsive neurons
subplot(2,2,1);
pie([group_summary.stim1_responsive_fdr, group_summary.stim2_responsive_fdr, ...
    group_summary.total_neurons - (group_summary.stim1_responsive_fdr + group_summary.stim2_responsive_fdr)]);
legend({'Safe Stimulus Preferring', 'Different Stimulus Preferring', 'Non-selective'}, 'Location', 'southoutside');
title('Proportion of Responsive Neurons (FDR-corrected)');

% Plot 2: Average calcium event rates
subplot(2,2,2);
bar([mean(all_mean_stim1), mean(all_mean_stim2)]);
hold on;
errorbar([1, 2], [mean(all_mean_stim1), mean(all_mean_stim2)], ...
         [std(all_mean_stim1)/sqrt(length(all_mean_stim1)), std(all_mean_stim2)/sqrt(length(all_mean_stim2))], ...
         'k', 'LineStyle', 'none');
hold off;
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Safe Stimulus', 'Different Stimulus'});
ylabel('Mean Event Rate (events/sec)');
title('Average Calcium Event Rate');

% Plot 3: Scatter plot of event rates for all neurons
subplot(2,2,[3,4]);
scatter(all_mean_stim1, all_mean_stim2, 25, all_responsive_fdr, 'filled');
hold on;
plot([0, max([all_mean_stim1; all_mean_stim2])], [0, max([all_mean_stim1; all_mean_stim2])], 'k--');
hold off;
colormap([0 0 1; 1 0 0; 0.7 0.7 0.7]);  % Blue for stim1, Red for stim2, Gray for non-selective
xlabel('Safe Stimulus Event Rate (events/sec)');
ylabel('Different Stimulus Event Rate (events/sec)');
title('Stimulus Preference by Neuron');

% Save the results
save('calcium_events_analysis.mat', 'calcium_events', 'group_summary');