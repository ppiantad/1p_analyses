%data are from pan-neuronal BLA imaging and SLEAP. Load these final .mat
%files as necessary!


%%
% uncomment below to load saved data if starting fresh (only uncomment
% whichever are necessary):

% load('BLA_panneuronal_Risk_2024_01_04.mat')
% load(BLA_panneuronal_SLEAP_data_2024_01_26.mat')
% load('events_panneuronal_RDT_Large_Pre-Choice+Large_Post-Choice_No_SHK+SHK.mat')

action_activated_alone = respClass_all_array{1, 1} == 1 & respClass_all_array{1, 2} == 3 & respClass_all_array{1, 3} == 3
sum(action_activated_alone)

consumption_activated_alone = respClass_all_array{1, 1} == 3 & respClass_all_array{1, 2} == 1 & respClass_all_array{1, 3} == 3
sum(consumption_activated_alone)

Shock_activated_alone = respClass_all_array{1, 1} == 3 & respClass_all_array{1, 2} == 3 & respClass_all_array{1, 3} == 1
sum(Shock_activated_alone)


%% Run eventRelatedActivity first for whatever events you want to identify responsive neurons for

clear zall_first_event zall_second_event zall_third_event


select_mouse = 'BLA_Insc_32';

select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'RDT_D1';

second_session = 'RDT_D1';


figure; plot(ts1, mean(neuron_mean_mouse{select_mouse_index, 1}(respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_Minus_4to0.REW_Large_BLOCK_1 == 1, :)))
hold on; plot(ts1, mean(neuron_mean_mouse{select_mouse_index, 2}(respClass_mouse.(select_mouse).(first_session).collectionTime.Outcome_1to3.REW_Large_SHK_0  == 1, :)))
hold on; plot(ts1, mean(neuron_mean_mouse{select_mouse_index, 3}(respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_0to2.SHK_1 == 1, :)))
% hold on; plot(ts1, mean(final_SLEAP.(select_mouse).(first_session).choiceTime.unitXTrials.zall_motion(final_SLEAP.(select_mouse).(first_session).BehavData.bigSmall == 1.2, :)))


action_and_shock = sum(respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_Minus_4to0.REW_Large_BLOCK_1 == 1 & respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_0to2.SHK_1 == 1)
consumption_and_shock = sum(respClass_mouse.(select_mouse).(first_session).collectionTime.Outcome_1to3.REW_Large_SHK_0 == 1 & respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_0to2.SHK_1 == 1)
action_and_consumption = sum(respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_Minus_4to0.REW_Large_BLOCK_1 == 1 & respClass_mouse.(select_mouse).(first_session).collectionTime.Outcome_1to3.REW_Large_SHK_0 == 1)

action_alone = respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_Minus_4to0.REW_Large_BLOCK_1 == 1 & respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_0to2.SHK_1 == 3 & respClass_mouse.(select_mouse).(first_session).collectionTime.Outcome_1to3.REW_Large_SHK_0 == 3
sum(action_alone)

consumption_alone = respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_Minus_4to0.REW_Large_BLOCK_1 == 3 & respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_0to2.SHK_1 == 3 & respClass_mouse.(select_mouse).(first_session).collectionTime.Outcome_1to3.REW_Large_SHK_0 == 1
sum(consumption_alone)

shock_alone = respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_Minus_4to0.REW_Large_BLOCK_1 == 3 & respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_0to2.SHK_1 == 1 & respClass_mouse.(select_mouse).(first_session).collectionTime.Outcome_1to3.REW_Large_SHK_0 == 3
sum(consumption_alone)


%% organize data - lots of this is manual due to need to specify which events have been previously filtered


caTraceTrials = final.(select_mouse).(first_session).CNMFe_data.C_raw(respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_Minus_4to0.REW_Large_BLOCK_1 == 1, :);

for h = 1:size(caTraceTrials,1)
    zb(h) = mean(caTraceTrials(h,:)); %baseline mean
    zsd(h) = std(caTraceTrials(h,:)); %baseline std
    tmp = 0;
    for j = 1:size(caTraceTrials,2)
        tmp = tmp+1;
        zall(h,tmp) = (caTraceTrials(h,tmp) - zb(h))/zsd(h);
        % Bound the normalized values between -1 and 1
        % zall(h,tmp) = max(-1, min(1, zall(h,tmp)));

    end

end
for z = 1:size(zall, 1)
    % Apply Savitzky-Golay filter to each row
    zall(z, :) = sgolayfilt(zall(z, :), 9, 21);
end
zall_first_event = zall; 
clear zall caTraceTrials

caTraceTrials = final.(select_mouse).(first_session).CNMFe_data.C_raw(respClass_mouse.(select_mouse).(first_session).collectionTime.Outcome_1to3.REW_Large_SHK_0 == 1, :);

for h = 1:size(caTraceTrials,1)
    zb(h) = mean(caTraceTrials(h,:)); %baseline mean
    zsd(h) = std(caTraceTrials(h,:)); %baseline std
    tmp = 0;
    for j = 1:size(caTraceTrials,2)
        tmp = tmp+1;
        zall(h,tmp) = (caTraceTrials(h,tmp) - zb(h))/zsd(h);
        % Bound the normalized values between -1 and 1
        % zall(h,tmp) = max(-1, min(1, zall(h,tmp)));
    end

end
for z = 1:size(zall, 1)
    % Apply Savitzky-Golay filter to each row
    zall(z, :) = sgolayfilt(zall(z, :), 9, 21);
end
zall_second_event = zall;
clear zall caTraceTrials


caTraceTrials = final.(select_mouse).(first_session).CNMFe_data.C_raw(respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_0to2.SHK_1  == 1, :);

for h = 1:size(caTraceTrials,1)
    zb(h) = mean(caTraceTrials(h,:)); %baseline mean
    zsd(h) = std(caTraceTrials(h,:)); %baseline std
    tmp = 0;
    for j = 1:size(caTraceTrials,2)
        tmp = tmp+1;
        zall(h,tmp) = (caTraceTrials(h,tmp) - zb(h))/zsd(h);
        % Bound the normalized values between -1 and 1
        % zall(h,tmp) = max(-1, min(1, zall(h,tmp)));

    end

end
for z = 1:size(zall, 1)
    % Apply Savitzky-Golay filter to each row
    zall(z, :) = sgolayfilt(zall(z, :), 9, 21);
end
zall_third_event = zall; 
clear zall caTraceTrials




caTraceTrials = final_SLEAP.(select_mouse).(first_session).SLEAP_data.vel_cm_s';


zb= mean(caTraceTrials); %baseline mean
zsd = std(caTraceTrials); %baseline std
tmp = 0;
for j = 1:size(caTraceTrials,2)
    tmp = tmp+1;
    zall(tmp) = (caTraceTrials(tmp) - zb)/zsd;
    % Bound the normalized values between -1 and 1
    % zall(h,tmp) = max(-1, min(1, zall(h,tmp)));
end

for z = 1:size(zall, 1)
    % Apply Savitzky-Golay filter to each row
    zall(z, :) = sgolayfilt(zall(z, :), 9, 21);
end
zall_motion = zall;
clear zall caTraceTrials
%%
% in order to trim off excess (because calcium recording starts before
% behavior), you'll need the start time. unfortunately I don't save that
% variable in the "final" struct, but you can get it from the adjustment to
% some of the columns of BehavData. e.g., see below - but make sure to
% update the session etc as necessary! 

BehavData = final.(select_mouse).(first_session).choiceTime.uv.BehavData;
% because the first trial possible is ALWAYS 60 seconds after ABET is
% issued, you can determine what adjustment has been made to this column
% (adding time to account for calcium recording starting first) by
% subtracting off 60 from the first element
stTime = BehavData.TrialPossible(1)-60; 

%%
num_samples = size(zall_second_event, 2)
sampling_frequency = (final.(select_mouse).(first_session).choiceTime.uv.dt)*100;
% Create a time array
time_array = (0:(num_samples-1)) / sampling_frequency;

time_columns_to_add = (0:(num_samples-1)) / sampling_frequency;


trim_length = size(final_SLEAP.(select_mouse).(first_session).zscored_SLEAP_data_velocity, 2);

time_array = time_array(:, 1:trim_length);

% it is necessary to add some empty columns because the way it is right
% now, SLEAP data starts when ABET starts. might make more sense to just
% trim off the initial samples of the session-long calcium data instead
% though? 
num_columns_added_to_SLEAP = round(sampling_frequency * stTime);

buffer_to_SLEAP = zeros(1, num_columns_added_to_SLEAP);

zall_motion = [buffer_to_SLEAP zall_motion]




figure; plot(time_array, mean(zall_first_event(:, 1:trim_length)))
hold on; plot(time_array, mean(zall_second_event(:, 1:trim_length)))
hold on; plot(time_array, mean(zall_third_event(:, 1:trim_length)))

% hold on; plot(time_array, zall_motion(:, 1:trim_length))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
xline(BehavData.choiceTime(BehavData.shock == 1), '--y')

%%
% Assuming zall_first_event and zall_second_event are your 22x19903 arrays
[row_count_first, column_count_first] = size(zall_first_event);
[row_count_second, column_count_second] = size(zall_second_event);

% Select neurons from each array
neurons_to_plot = 3;
selected_neurons_first = zall_first_event(1:neurons_to_plot, 1:trim_length);
selected_neurons_second = zall_second_event(1:neurons_to_plot, 1:trim_length);
selected_neurons_third = zall_third_event(1:neurons_to_plot, 1:trim_length);

% Combine selected neurons into one array
combined_neurons = [selected_neurons_first; selected_neurons_second; selected_neurons_third];



% Define spacing and line thickness
row_spacing = 5; % Adjust as needed
line_thickness = 1.0; % Adjust as needed

% Create a figure
figure;

% Loop through each row and plot the neural activity for the combined neurons
for i = 1:size(combined_neurons, 1)
    % Extract neural activity for the current neuron
    neuron_activity = combined_neurons(i, :);
    
    % Calculate the y-coordinate for the current neuron's plot
    y_coordinate = (i - 1) * row_spacing;
    
    % Plot the neural activity with the specified x and y coordinates
    plot(time_array, neuron_activity + y_coordinate, 'k', 'LineWidth', line_thickness);
    
    hold on; % To overlay all plots on the same figure
end

% Customize the plot
title('Neural Activity of Selected Neurons');
xlabel('Time'); % Change as per your time unit
ylabel('Neuron #'); % Customize y-axis title

grid off; % Add grid lines

% Remove ticks and labels on the top and right axes
set(gca, 'TickDir', 'out');
set(gca, 'Box', 'off');
ax = gca;

xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')

xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
xline(BehavData.choiceTime(BehavData.shock == 1), '--y')

% Show legend or additional information (if needed)

% Hold off to stop overlaying subsequent plots on the same figure
hold off;