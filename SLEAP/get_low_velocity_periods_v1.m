
uv.evtWin = [-8 8]; %what time do you want to look at around each event [-2 8] [-10 5] [-10 10]
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
session_to_analyze = 'RDT_D1';


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





%%

animalIDs = (fieldnames(final));
block_categorized_data = struct(); % Create structure to store block-categorized data

for ii = 1:size(animalIDs,1)
    currentanimal = char(animalIDs(ii));
    
    if isfield(final.(currentanimal), session_to_analyze)
        BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;
        
        % Step 1: Create a logical index for velocity < 1
        velocity_threshold = 1;
        velocity_below_threshold = final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data.vel_cm_s < velocity_threshold;
        
        % Step 2: Identify start indices of low-velocity bouts
        start_indices = find([0; diff(velocity_below_threshold)] == 1);
        
        % Step 3: Filter out bouts where the velocity stays < 1 for at least one second,
        % the separation between bouts is > 2 seconds, and not within ±1 second of collectionTime
        min_duration_seconds = 1;
        min_separation_seconds = 1;
        min_mean_velocity_before_start = 2;
        valid_start_indices = [];
        
        for i = 1:length(start_indices)
            current_start_index = start_indices(i);
            current_end_index = min(length(velocity_below_threshold), current_start_index + 10 * min_duration_seconds);
            
            if all(velocity_below_threshold(current_start_index:current_end_index))
                % Check if the separation between bouts is > 2 seconds
                if i == 1 || (final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data.idx_time(start_indices(i)) - final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data.idx_time(start_indices(i-1)) > min_separation_seconds)
                    % Check if not within ±1 second of collectionTime
                    current_start_time = final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data.idx_time(current_start_index);
                    collection_times = final_SLEAP.(currentanimal).(session_to_analyze).BehavData.collectionTime;
                    
                    if ~any(abs(current_start_time - collection_times) <= 5)
                        % Check if mean velocity in the 1 second before start index is > 2 cm/s
                        if current_start_index > 10 % Make sure we have enough previous data points
                            mean_velocity_before_start = mean(final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data.vel_cm_s(current_start_index - 10:current_start_index - 1));
                            if mean_velocity_before_start > min_mean_velocity_before_start
                                valid_start_indices = [valid_start_indices, current_start_index];
                            end
                        end
                    end
                end
            end
        end
        
        % Step 4: Extract corresponding values from 'idx_time'
        periods_with_low_velocity{ii} = final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data.idx_time(valid_start_indices);
        
        % NEW CODE: Determine Block start and end times
        blocks = unique(BehavData.Block);
        block_times = struct();
        
        % Initialize arrays to store categorized periods
        block_categorized_data.(currentanimal).block1_periods = [];
        block_categorized_data.(currentanimal).block2_periods = [];
        block_categorized_data.(currentanimal).block3_periods = [];
        block_categorized_data.(currentanimal).outside_blocks = [];
        
        % Get the start and end times for each block
        for b = 1:length(blocks)
            current_block = blocks(b);
            
            % Find indices where Block equals current_block
            block_indices = find(BehavData.Block == current_block);
            
            if ~isempty(block_indices)
                % Get start time from the first trial of the block
                block_start_time = BehavData.stTime(block_indices(1));
                
                % Get end time from the last trial of the block
                block_end_time = BehavData.collectionTime(block_indices(end));
                
                % Store the times
                block_times.(['block', num2str(current_block)]) = [block_start_time, block_end_time];
                
                fprintf('Animal %s - Block %d: Start=%.2f, End=%.2f\n', currentanimal, current_block, block_start_time, block_end_time);
            end
        end
        
        % Categorize periods of low velocity based on which block they fall into
        if ~isempty(periods_with_low_velocity{ii})
            for p = 1:length(periods_with_low_velocity{ii})
                current_period = periods_with_low_velocity{ii}(p);
                categorized = false;
                
                % Check which block this period belongs to
                for b = 1:length(blocks)
                    current_block = blocks(b);
                    block_key = ['block', num2str(current_block)];
                    
                    if isfield(block_times, block_key)
                        block_start = block_times.(block_key)(1);
                        block_end = block_times.(block_key)(2);
                        
                        % If the period is within this block's time range
                        if current_period >= block_start && current_period <= block_end
                            block_categorized_data.(currentanimal).([block_key, '_periods']) = [block_categorized_data.(currentanimal).([block_key, '_periods']), current_period];
                            categorized = true;
                            break;
                        end
                    end
                end
                
                % If not within any block, categorize as outside blocks
                if ~categorized
                    block_categorized_data.(currentanimal).outside_blocks = [block_categorized_data.(currentanimal).outside_blocks, current_period];
                end
            end
        end
        
        % Display summary
        fprintf('\nAnimal %s - Summary of low velocity periods:\n', currentanimal);
        fprintf('Block 1: %d periods\n', length(block_categorized_data.(currentanimal).block1_periods));
        fprintf('Block 2: %d periods\n', length(block_categorized_data.(currentanimal).block2_periods));
        fprintf('Block 3: %d periods\n', length(block_categorized_data.(currentanimal).block3_periods));
        fprintf('Outside blocks: %d periods\n', length(block_categorized_data.(currentanimal).outside_blocks));
        
        % You could also compute statistics, like density of periods per block
        if isfield(block_times, 'block1') && ~isempty(block_times.block1)
            block1_duration = block_times.block1(2) - block_times.block1(1);
            block1_density = length(block_categorized_data.(currentanimal).block1_periods) / block1_duration;
            fprintf('Block 1 density: %.4f periods/second\n', block1_density);
        end
        
        if isfield(block_times, 'block2') && ~isempty(block_times.block2)
            block2_duration = block_times.block2(2) - block_times.block2(1);
            block2_density = length(block_categorized_data.(currentanimal).block2_periods) / block2_duration;
            fprintf('Block 2 density: %.4f periods/second\n', block2_density);
        end
        
        if isfield(block_times, 'block3') && ~isempty(block_times.block3)
            block3_duration = block_times.block3(2) - block_times.block3(1);
            block3_density = length(block_categorized_data.(currentanimal).block3_periods) / block3_duration;
            fprintf('Block 3 density: %.4f periods/second\n', block3_density);
        end
    end
end

% Now block_categorized_data contains all the low velocity periods organized by animal and trial block
disp('Analysis complete!')

%%
% After your previous code completes, add this to create the bar plot

% Prepare data for plotting
n_animals = size(animalIDs, 1);
block1_counts = zeros(n_animals, 1);
block2_counts = zeros(n_animals, 1);
block3_counts = zeros(n_animals, 1);

% Extract counts for each animal and block
for ii = 1:n_animals
    currentanimal = char(animalIDs(ii));
    
    % Check if this animal has data
    if isfield(block_categorized_data, currentanimal)
        if isfield(block_categorized_data.(currentanimal), 'block1_periods')
            block1_counts(ii) = length(block_categorized_data.(currentanimal).block1_periods);
        end
        
        if isfield(block_categorized_data.(currentanimal), 'block2_periods')
            block2_counts(ii) = length(block_categorized_data.(currentanimal).block2_periods);
        end
        
        if isfield(block_categorized_data.(currentanimal), 'block3_periods')
            block3_counts(ii) = length(block_categorized_data.(currentanimal).block3_periods);
        end
    end
end

% Calculate means and SEMs
mean_counts = [mean(block1_counts), mean(block2_counts), mean(block3_counts)];
sem_counts = [std(block1_counts)/sqrt(n_animals), std(block2_counts)/sqrt(n_animals), std(block3_counts)/sqrt(n_animals)];

% Create figure
figure('Position', [100, 100, 800, 600]);
bar_h = bar(1:3, mean_counts, 0.6);
hold on;

% Set colors for the bars (you can customize)
bar_h.FaceColor = 'flat';
bar_h.CData(1,:) = [0.2, 0.6, 0.8]; % Block 1 color
bar_h.CData(2,:) = [0.8, 0.4, 0.2]; % Block 2 color
bar_h.CData(3,:) = [0.3, 0.7, 0.3]; % Block 3 color

% Add error bars
errorbar(1:3, mean_counts, sem_counts, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Add individual data points (scatter)
% Calculate x-coordinates for scatter points with some jitter
scatter_x1 = ones(size(block1_counts)) + 0.1*(rand(size(block1_counts))-0.5);
scatter_x2 = 2*ones(size(block2_counts)) + 0.1*(rand(size(block2_counts))-0.5);
scatter_x3 = 3*ones(size(block3_counts)) + 0.1*(rand(size(block3_counts))-0.5);

% Plot scatter points
scatter(scatter_x1, block1_counts, 70, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
scatter(scatter_x2, block2_counts, 70, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
scatter(scatter_x3, block3_counts, 70, 'k', 'filled', 'MarkerFaceAlpha', 0.7);

% Add lines connecting the same animal across blocks
for ii = 1:n_animals
    plot([scatter_x1(ii), scatter_x2(ii), scatter_x3(ii)], ...
         [block1_counts(ii), block2_counts(ii), block3_counts(ii)], ...
         'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
end

% Customize the plot
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'Box', 'off');
set(gca, 'XTick', 1:3, 'XTickLabel', {'Block 1', 'Block 2', 'Block 3'});
ylabel('Sum low velocity periods', 'FontSize', 16);
% title('Mean Low Velocity Periods by Trial Block', 'FontSize', 18);

% Add significance test information
[h_12, p_12] = ttest(block1_counts, block2_counts);
[h_13, p_13] = ttest(block1_counts, block3_counts);
[h_23, p_23] = ttest(block2_counts, block3_counts);

% Display significance test results in command window
fprintf('\nStatistical Tests (paired t-test):\n');
fprintf('Block 1 vs Block 2: p = %.4f\n', p_12);
fprintf('Block 1 vs Block 3: p = %.4f\n', p_13);
fprintf('Block 2 vs Block 3: p = %.4f\n', p_23);

% Add a legend for clarity
% legend(bar_h, {'Block 1', 'Block 2', 'Block 3'}, 'Location', 'best', 'FontSize', 12);

% Optional: Add statistical significance markers on the plot
y_max = max(mean_counts + sem_counts) * 1.1;
if p_12 < 0.05
    line([1, 2], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.5);
    if p_12 < 0.001
        text(1.5, y_max*1.02, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p_12 < 0.01
        text(1.5, y_max*1.02, '**', 'HorizontalAlignment', 'center', 'FontSize', 14);
    else
        text(1.5, y_max*1.02, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    end
    y_max = y_max * 1.05;
end

if p_13 < 0.05
    line([1, 3], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.5);
    if p_13 < 0.001
        text(2, y_max*1.02, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p_13 < 0.01
        text(2, y_max*1.02, '**', 'HorizontalAlignment', 'center', 'FontSize', 14);
    else
        text(2, y_max*1.02, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    end
    y_max = y_max * 1.05;
end

if p_23 < 0.05
    line([2, 3], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.5);
    if p_23 < 0.001
        text(2.5, y_max*1.02, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p_23 < 0.01
        text(2.5, y_max*1.02, '**', 'HorizontalAlignment', 'center', 'FontSize', 14);
    else
        text(2.5, y_max*1.02, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    end
end

% % Save the figure if desired
% saveas(gcf, 'low_velocity_periods_by_block.png');
% saveas(gcf, 'low_velocity_periods_by_block.fig');

% Optional: Create a table with the data for export
results_table = table(animalIDs, block1_counts, block2_counts, block3_counts, ...
    'VariableNames', {'AnimalID', 'Block1_Counts', 'Block2_Counts', 'Block3_Counts'});
% writetable(results_table, 'low_velocity_periods_by_block.csv');

%%


for i = 1 %could loop through multiple mice like this if you had it
    eTS = periods_with_low_velocity; %BehavData.(uv.behav); %get time stamps
    velocity = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data.vel_cm_s; %get velocity
    %     ca = neuron.S; %get binarized calcium
    
    velocity_time = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data.idx_time'; % time trace
    % velocity_time  = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data.idx_time(1):uv.dt:length(velocity)*uv.dt; %generate time trace

    %calculate time windows for each event
    evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
    numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate
    %
    tic

    % initialize trial matrices
    velocity_TraceTrials = NaN(size(eTS,1),numMeasurements); %
    % unitTrace = velocity(u,:); %get trace
    %             %%
    for t = 1:size(eTS,1)
        % set each trial's temporal boundaries
        timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
        BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
        if min(timeWin) > min(velocity_time ) & max(timeWin) < max(velocity_time )    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
            % get unit event counts in trials
            % get unit ca traces in trials
            idx = velocity_time > min(timeWin) & velocity_time <= max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
            sum(idx)
            bl_idx = velocity_time > min(BL_win) & velocity_time <= max(BL_win);
            %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
            velocity_TraceTrials(t,1:sum(idx)) = velocity(idx);
            % zb(t,:) = mean(unitTrace(bl_idx)); %baseline mean
            % zsd(t,:) = std(unitTrace(bl_idx)); %baseline std
            velocity_zb(t,:) = mean(velocity_TraceTrials(t,:)); %baseline mean
            velocity_zsd(t,:) = std(velocity_TraceTrials(t,:)); %baseline std


            tmp = 0;
            for j = 1:size(velocity_TraceTrials,2)
                tmp = tmp+1;
                velocity_zall(t,tmp) = (velocity_TraceTrials(t,j) - velocity_zb(t))/velocity_zsd(t);
            end
            clear j;
        elseif ~(min(timeWin) > min(velocity_time ) & max(timeWin) < max(velocity_time ))
            continue


        end
        clear idx timeWin BL_win bl_idx

        unitXTrials.velocity_Traces = velocity_TraceTrials;
        unitXTrials.velocity_zb = velocity_zb;
        unitXTrials.velocity_zsd = velocity_zsd;
        % unitXTrials(u).zall = velocity_zall;

        % store unit averaged data
        % unitAVG.caTraces(u,:) = nanmean(velocity_TraceTrials);           %store trial averaged calcium traces
        % unitSEM.caTraces(u,:) = std(velocity_TraceTrials,'omitnan')/sqrt(size(velocity_TraceTrials,1));
        clear caEvtCtTrials caEvtRateTrials unitTrace idx
    end
end