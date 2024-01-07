%% Run eventRelatedActivity first for whatever events you want to identify responsive neurons for

select_mouse = 'BLA_Insc_27';

select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'RM_D1';

second_session = 'RM_D1';

action_activated_data = neuron_mean_mouse_unnormalized{select_mouse_index, 2}(respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_Minus_4to0.REW_Large.activated== 1, :);
consumption_activated_data = neuron_mean_mouse_unnormalized{select_mouse_index, 3}(respClass_mouse.(select_mouse).(first_session).collectionTime.Outcome_1to3.REW_Large.activated == 1, :);
velocity_data = final_SLEAP.(select_mouse).(first_session).choiceTime.unitXTrials.velocity_trace_trials(final_SLEAP.(select_mouse).(first_session).BehavData.bigSmall == 1.2, :);


test_data = [action_activated_data(:, 1:end-1); consumption_activated_data(:, 1:end-1); velocity_data];

zscore_data = zscore(test_data, 0 , 2);




figure; plot(ts1, mean(neuron_mean_mouse{select_mouse_index, 2}(respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_Minus_4to0.REW_Large.activated== 1, :)))
hold on; plot(ts1, mean(neuron_mean_mouse{select_mouse_index, 3}(respClass_mouse.(select_mouse).(first_session).collectionTime.Outcome_1to3.REW_Large.activated == 1, :)))
hold on; plot(ts1, mean(final_SLEAP.(select_mouse).(first_session).choiceTime.unitXTrials.zall_motion(final_SLEAP.(select_mouse).(first_session).BehavData.bigSmall == 1.2, :)))


%% organize data - lots of this is manual due to need to specify which events have been previously filtered


caTraceTrials = final.(select_mouse).(first_session).CNMFe_data.C_raw(respClass_mouse.(select_mouse).(first_session).choiceTime.Outcome_Minus_4to0.REW_Large.activated== 1, :);

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

caTraceTrials = final.(select_mouse).(first_session).CNMFe_data.C_raw(respClass_mouse.(select_mouse).(first_session).collectionTime.Outcome_1to3.REW_Large.activated== 1, :);

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
hold on; plot(time_array, zall_motion(:, 1:trim_length))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2), '--k')




