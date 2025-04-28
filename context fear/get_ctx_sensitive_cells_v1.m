
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

ca_data_type = "S"; % C % C_raw %S
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes
% CNMFe_data.spike_prob: CASCADE inferred spikes - multiply x sampling rate
% (10) for spike rate


session_to_analyze = 'D4';

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
neuron_num_experimental = 0;
neuron_num_one_context = 0;
neuron_num_no_shock = 0;

animalIDs = (fieldnames(final));


% Initialize a structure to store calcium event counts
calcium_events = struct();

% Process each animal
animalIDs = fieldnames(final);
for ii = 1:size(animalIDs,1)
    currentanimal = char(animalIDs(ii));
    
    % Check if animal has the session we want to analyze
    if isfield(final.(currentanimal), session_to_analyze)
        current_animal_treatment{ii} = final.(currentanimal).experimental_grp;
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
        calcium_events.(currentanimal).treatment_condition = current_animal_treatment{ii};
        
        % Process each neuron
        for neuron_idx = 1:num_neurons
            neuron_num = neuron_num + 1;
            if strcmp(current_animal_treatment{ii}, 'Experimental')
                neuron_num_experimental = neuron_num_experimental + 1;
            elseif strcmp(current_animal_treatment{ii}, 'One Context')
                neuron_num_one_context = neuron_num_one_context + 1;

            elseif strcmp(current_animal_treatment{ii}, 'No Shock')
                neuron_num_no_shock = neuron_num_no_shock + 1;

            end

            neuron_trace = ca(neuron_idx, :);
            
            % Find peaks in the trace (calcium events)
            % [peaks, peak_locs] = findpeaks(neuron_trace, 'MinPeakHeight', mean(neuron_trace) + 1*std(neuron_trace), ...
            %                               'MinPeakDistance', 3);


            [peaks, peak_locs] = findpeaks(neuron_trace);

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

            % norm_stim1_events = calcium_events.(currentanimal).stimulus_type_1(neuron_idx, :);
            % norm_stim2_events = calcium_events.(currentanimal).stimulus_type_2(neuron_idx, :);
            
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
% [fdr_p, fdr_mask] = mafdr(all_p_values, 'BHFDR', false);
% all_responsive_fdr = zeros(size(all_responsive));
% for i = 1:length(fdr_mask)
%     if fdr_mask(i)
%         if all_mean_stim1(i) > all_mean_stim2(i)
%             all_responsive_fdr(i) = 1;
%         else
%             all_responsive_fdr(i) = 2;
%         end
%     end
% end

% Apply multiple comparisons correction (False Discovery Rate)
[fdr_p, fdr_mask] = mafdr(all_p_values, 'BHFDR', false);
all_responsive_fdr = zeros(size(all_responsive));
for i = 1:length(fdr_mask)
    fdr_mask(i);
    if all_mean_stim1(i) > all_mean_stim2(i)
        all_responsive_fdr(i) = 1;
    else
        all_responsive_fdr(i) = 2;
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


experimental_mice_num = 0;
one_context_mice_num = 0;
no_shock_mice_num = 0;

for ii = 1:length(animalIDs)
    currentanimal = char(animalIDs(ii));
    currentTreatment = calcium_events.(currentanimal).treatment_condition;
    if strcmp(currentTreatment, 'Experimental')
        experimental_mice_num = experimental_mice_num + 1
        percent_responsive_experimental_mice(experimental_mice_num) = calcium_events.(currentanimal).percent_responsive;
        percent_increase_first_third_fifth_experimental_mice(experimental_mice_num) = calcium_events.(currentanimal).percent_stim1_responsive;
        percent_increase_second_fourth_sixth_experimental_mice(experimental_mice_num) = calcium_events.(currentanimal).percent_stim2_responsive;
    elseif strcmp(currentTreatment, 'One Context')
        one_context_mice_num = one_context_mice_num + 1;
        percent_responsive_one_context_mice(one_context_mice_num) = calcium_events.(currentanimal).percent_responsive;
        percent_increase_first_third_fifth_context_mice(one_context_mice_num) = calcium_events.(currentanimal).percent_stim1_responsive;
        percent_increase_second_fourth_sixth_context_mice(one_context_mice_num) = calcium_events.(currentanimal).percent_stim2_responsive;
    elseif strcmp(currentTreatment, 'No Shock')
        no_shock_mice_num = no_shock_mice_num + 1;
        percent_responsive_no_shock_mice(no_shock_mice_num) = calcium_events.(currentanimal).percent_responsive;
        percent_increase_first_third_fifth_no_shock_mice(no_shock_mice_num) = calcium_events.(currentanimal).percent_stim1_responsive;
        percent_increase_second_fourth_sixth_no_shock_mice(no_shock_mice_num) = calcium_events.(currentanimal).percent_stim2_responsive;
    end
end

% Calculate means
means = [mean(percent_responsive_experimental_mice), ...
         mean(percent_responsive_one_context_mice), ...
         mean(percent_responsive_no_shock_mice)];

% Calculate SEM (Standard Error of the Mean)
sem = [std(percent_responsive_experimental_mice)/sqrt(length(percent_responsive_experimental_mice)), ...
       std(percent_responsive_one_context_mice)/sqrt(length(percent_responsive_one_context_mice)), ...
       std(percent_responsive_no_shock_mice)/sqrt(length(percent_responsive_no_shock_mice))];

% Create a figure
figure('Position', [100, 100, 600, 500]);

% Create bar plot
b = bar(means, 'FaceColor', 'flat');
b.CData(1,:) = [0.3, 0.6, 0.8]; % Color for experimental
b.CData(2,:) = [0.8, 0.4, 0.3]; % Color for one context
b.CData(3,:) = [0.4, 0.7, 0.4]; % Color for no shock

hold on;

% Add error bars
errorbar(1:3, means, sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Add individual data points
% Need to calculate x-positions for the scatter points (with some jitter)
x1 = repmat(1, size(percent_responsive_experimental_mice)) + (rand(size(percent_responsive_experimental_mice))-0.5)*0.2;
x2 = repmat(2, size(percent_responsive_one_context_mice)) + (rand(size(percent_responsive_one_context_mice))-0.5)*0.2;
x3 = repmat(3, size(percent_responsive_no_shock_mice)) + (rand(size(percent_responsive_no_shock_mice))-0.5)*0.2;

% Plot scatter points
scatter(x1, percent_responsive_experimental_mice, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
scatter(x2, percent_responsive_one_context_mice, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
scatter(x3, percent_responsive_no_shock_mice, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.7);

% Set axis labels and title
% xlabel('Condition', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Percent responsive (%)', 'FontSize', 14, 'FontWeight', 'bold');
% title('Percent Responsive Mice by Condition', 'FontSize', 16, 'FontWeight', 'bold');

% Set x-tick labels
xticks(1:3);
xticklabels({'Experimental', 'One Context', 'No Shock'});
xtickangle(45);

% Improve appearance
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'off');
ylim([0, max(means) + max(sem) + 0.1]);

% Add a legend (if needed)
legend('Experimental', 'One Context', 'No Shock', 'Location', 'NorthEast');

% Adjust spacing
% grid on;
% grid minor;
hold off;


%%
% After completing the analysis above, create a combined visualization of neuron activity
figure('Name', 'Significant and Highly Active Neurons', 'Position', [100, 100, 1400, 900]);

% Create time vector for x-axis (in seconds)
time_vector = (0:(total_time_points-1)) / sampling_rate;

% Initialize arrays to store neuron information with their p-values and mean activity
stim1_neurons_data = struct('animal', {}, 'neuron_idx', {}, 'p_value', {}, 'mean_activity', {}, 'trace', {});
stim2_neurons_data = struct('animal', {}, 'neuron_idx', {}, 'p_value', {}, 'mean_activity', {}, 'trace', {});

% Define significance and activity thresholds
p_value_threshold = 0.05;
stim1_activity_threshold = 0.10;  % For safe-active neurons
stim2_activity_threshold = 0.10;  % For different-active neurons

for ii = 1:size(animalIDs,1)
    currentanimal = char(animalIDs(ii));
    
    % Check if animal has the session we want to analyze
    if isfield(final.(currentanimal), session_to_analyze)
        % Get calcium data
        ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);
        if strcmp(ca_data_type, 'S')
            ca = full(ca);
        end
        
        % Get responsive neurons, p-values, and mean activity during stimuli
        responsive_neurons = calcium_events.(currentanimal).responsive_neurons;
        neuron_p_values = calcium_events.(currentanimal).p_values;
        mean_events_stim1 = calcium_events.(currentanimal).mean_events_stim1;  % Assuming this exists
        mean_events_stim2 = calcium_events.(currentanimal).mean_events_stim2;  % Assuming this exists
        
        % Extract stim1 preferring neurons (safe-active)
        stim1_neurons = find(responsive_neurons == 1);
        for n_idx = 1:length(stim1_neurons)
            neuron_idx = stim1_neurons(n_idx);
            
            % Check if this neuron meets both criteria
            if neuron_p_values(neuron_idx) < p_value_threshold && mean_events_stim1(neuron_idx) > stim1_activity_threshold
                % Ensure trace is the right length by trimming if necessary
                trace = ca(neuron_idx, 1:min(size(ca, 2), length(time_vector)));
                
                % Normalize trace to [0,1]
                normalized_trace = (trace - min(trace)) / (max(trace) - min(trace));
                
                % Store neuron data
                new_data = struct('animal', currentanimal, ...
                                 'neuron_idx', neuron_idx, ...
                                 'p_value', neuron_p_values(neuron_idx), ...
                                 'mean_activity', mean_events_stim1(neuron_idx), ...
                                 'trace', normalized_trace);
                stim1_neurons_data = [stim1_neurons_data; new_data];
            end
        end
        
        % Extract stim2 preferring neurons (different-active)
        stim2_neurons = find(responsive_neurons == 2);
        for n_idx = 1:length(stim2_neurons)
            neuron_idx = stim2_neurons(n_idx);
            
            % Check if this neuron meets both criteria
            if neuron_p_values(neuron_idx) < p_value_threshold && mean_events_stim2(neuron_idx) > stim2_activity_threshold
                % Ensure trace is the right length by trimming if necessary
                trace = ca(neuron_idx, 1:min(size(ca, 2), length(time_vector)));
                
                % Normalize trace to [0,1]
                normalized_trace = (trace - min(trace)) / (max(trace) - min(trace));
                
                % Store neuron data
                new_data = struct('animal', currentanimal, ...
                                 'neuron_idx', neuron_idx, ...
                                 'p_value', neuron_p_values(neuron_idx), ...
                                 'mean_activity', mean_events_stim2(neuron_idx), ...
                                 'trace', normalized_trace);
                stim2_neurons_data = [stim2_neurons_data; new_data];
            end
        end
    end
end

% Sort neurons by p-value (ascending order)
if ~isempty(stim1_neurons_data)
    [~, idx] = sort([stim1_neurons_data.p_value]);
    stim1_neurons_data = stim1_neurons_data(idx);
end

if ~isempty(stim2_neurons_data)
    [~, idx] = sort([stim2_neurons_data.p_value]);
    stim2_neurons_data = stim2_neurons_data(idx);
end

% Take only top 5 neurons (or all if fewer than 5)
top_n = 5;
stim1_top_neurons = min(top_n, length(stim1_neurons_data));
stim2_top_neurons = min(top_n, length(stim2_neurons_data));

% Prepare data for plotting top safe-preferring neurons
stim1_traces = [];
stim1_neuron_ids = {};

for i = 1:stim1_top_neurons
    stim1_traces = [stim1_traces; stim1_neurons_data(i).trace];
    stim1_neuron_ids = [stim1_neuron_ids; {sprintf('%s-N%d (p=%.3f, act=%.2f)', ...
                        stim1_neurons_data(i).animal, ...
                        stim1_neurons_data(i).neuron_idx, ...
                        stim1_neurons_data(i).p_value, ...
                        stim1_neurons_data(i).mean_activity)}];
end

% Prepare data for plotting top different-preferring neurons
stim2_traces = [];
stim2_neuron_ids = {};

for i = 1:stim2_top_neurons
    stim2_traces = [stim2_traces; stim2_neurons_data(i).trace];
    stim2_neuron_ids = [stim2_neuron_ids; {sprintf('%s-N%d (p=%.3f, act=%.2f)', ...
                        stim2_neurons_data(i).animal, ...
                        stim2_neurons_data(i).neuron_idx, ...
                        stim2_neurons_data(i).p_value, ...
                        stim2_neurons_data(i).mean_activity)}];
end

% Trim time vector if needed to match trace length
if ~isempty(stim1_traces)
    time_vector = time_vector(1:size(stim1_traces, 2));
elseif ~isempty(stim2_traces)
    time_vector = time_vector(1:size(stim2_traces, 2));
end

% Plot both types of neurons on the same figure
figure(1);
hold on;

% Define spacing and layout parameters
spacing = 1.2;           % Spacing between traces
central_gap = 2;         % Gap between the two groups of neurons
total_neurons = stim1_top_neurons + stim2_top_neurons;

% Plot safe-preferring neurons (blue, bottom)
for i = 1:stim1_top_neurons
    plot(time_vector, stim1_traces(i,:) + (i-1)*spacing, 'b', 'LineWidth', 1.5);
end

% Plot different-preferring neurons (red, top)
for i = 1:stim2_top_neurons
    % Position them after safe-preferring neurons plus a gap
    plot(time_vector, stim2_traces(i,:) + (i-1)*spacing + stim1_top_neurons*spacing + central_gap, 'r', 'LineWidth', 1.5);
end

% Add shading for stimulus periods across entire plot height
y_max = (total_neurons-1) * spacing + central_gap + 1; 

% Shade stimulus type 1 periods (light blue)
for stim_idx = 1:size(stimulus_frames{1}, 1)
    frame_start = stimulus_frames{1}(stim_idx, 1);
    frame_end = stimulus_frames{1}(stim_idx, 2);
    x_start = frame_start / sampling_rate;
    x_end = frame_end / sampling_rate;
    
    % Make sure the stimulus periods are within our plot range
    if x_start <= time_vector(end) && x_end >= time_vector(1)
        % Adjust to be within range
        x_start = max(x_start, time_vector(1));
        x_end = min(x_end, time_vector(end));
        
        % Create shaded region
        rectangle('Position', [x_start, 0, x_end-x_start, y_max], ...
                 'FaceColor', [0.8, 0.9, 1, 0.2], 'EdgeColor', 'none');
    end
end

% Shade stimulus type 2 periods (light red)
for stim_idx = 1:size(stimulus_frames{2}, 1)
    frame_start = stimulus_frames{2}(stim_idx, 1);
    frame_end = stimulus_frames{2}(stim_idx, 2);
    x_start = frame_start / sampling_rate;
    x_end = frame_end / sampling_rate;
    
    % Make sure the stimulus periods are within our plot range
    if x_start <= time_vector(end) && x_end >= time_vector(1)
        % Adjust to be within range
        x_start = max(x_start, time_vector(1));
        x_end = min(x_end, time_vector(end));
        
        % Create shaded region
        rectangle('Position', [x_start, 0, x_end-x_start, y_max], ...
                 'FaceColor', [1, 0.8, 0.8, 0.2], 'EdgeColor', 'none');
    end
end

% Create combined y-tick positions and labels
ytick_positions = [];
ytick_labels = {};

% Add y-ticks for safe-preferring neurons
for i = 1:stim1_top_neurons
    ytick_positions = [ytick_positions, (i-1)*spacing];
    ytick_labels = [ytick_labels, stim1_neuron_ids(i)];
end

% Add y-ticks for different-preferring neurons
for i = 1:stim2_top_neurons
    ytick_positions = [ytick_positions, (i-1)*spacing + stim1_top_neurons*spacing + central_gap];
    ytick_labels = [ytick_labels, stim2_neuron_ids(i)];
end

% Set y-tick labels
yticks(ytick_positions);
yticklabels(ytick_labels);

% Add a horizontal line to separate the two groups
if stim1_top_neurons > 0 && stim2_top_neurons > 0
    separator_y = stim1_top_neurons * spacing + central_gap/2;
    line([time_vector(1), time_vector(end)], [separator_y, separator_y], 'Color', 'k', 'LineStyle', '--');
end

% Create custom legend
h_safe = plot(NaN, NaN, 'b', 'LineWidth', 1.5);
h_diff = plot(NaN, NaN, 'r', 'LineWidth', 1.5);
h_stim1 = rectangle('Position', [0, 0, 0, 0], 'FaceColor', [0.8, 0.9, 1, 0.2], 'EdgeColor', 'none');
h_stim2 = rectangle('Position', [0, 0, 0, 0], 'FaceColor', [1, 0.8, 0.8, 0.2], 'EdgeColor', 'none');

% Create legend
legend([h_safe, h_diff, h_stim1, h_stim2], ...
       {'Safe-Preferring Neurons', 'Different-Preferring Neurons', 'Safe Stimulus', 'Different Stimulus'}, ...
       'Location', 'northoutside', 'Orientation', 'horizontal');

% Add labels and title
xlabel('Time (seconds)');
title(['Top Neurons with p < ' num2str(p_value_threshold) ...
       ' and Activity > ' num2str(stim1_activity_threshold) ...
       ' (' num2str(stim1_top_neurons) ' Safe, ' num2str(stim2_top_neurons) ' Different)']);

% Add group labels to clarify the two sections
if stim1_top_neurons > 0
    text(time_vector(1) - 0.05*range(time_vector), (stim1_top_neurons-1)*spacing/2, 'Safe-Preferring', ...
         'Color', 'b', 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'Rotation', 90);
end

if stim2_top_neurons > 0
    text(time_vector(1) - 0.05*range(time_vector), (stim1_top_neurons)*spacing + central_gap + (stim2_top_neurons-1)*spacing/2, 'Different-Preferring', ...
         'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'Rotation', 90);
end

% Set y-limits with a bit of padding
ylim([-0.5, y_max + 0.5]);

hold off;

% Display information about the selection criteria
fprintf(['Neuron selection criteria:\n' ...
         '- Safe-preferring: p < %.3f and mean_events_stim1 > %.3f\n' ...
         '- Different-preferring: p < %.3f and mean_events_stim2 > %.3f\n' ...
         '- Found %d safe-preferring and %d different-preferring neurons meeting criteria\n'], ...
         p_value_threshold, stim1_activity_threshold, ...
         p_value_threshold, stim2_activity_threshold, ...
         length(stim1_neurons_data), length(stim2_neurons_data));

% Save figure
% savefig(1, 'significant_active_neurons.fig');
% saveas(1, 'significant_active_neurons.png');