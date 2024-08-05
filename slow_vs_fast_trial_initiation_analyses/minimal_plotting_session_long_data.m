animal_index = 7;
neuron_to_display = 72;
animalID_to_use = animalIDs{animal_index};
session_to_display = 'Pre_RDT_RM';

figure; imagesc(ts1, [], zall_mouse{animal_index, 3}{1, neuron_to_display})
figure; plot(ts1, mean(zall_mouse{animal_index, 3}{1, neuron_to_display}))


num_samples = size(final.(animalID_to_use).(session_to_display).CNMFe_data.C_raw, 2)
% sampling_frequency = (final.(select_mouse).(first_session).choiceTime.uv.dt)*100;
sampling_frequency = (final.(animalID_to_use).(session_to_display).uv.dt)*100;
% Create a time array
time_array = (0:(num_samples-1)) / sampling_frequency;

BehavData = final.(animalID_to_use).(session_to_display).uv.BehavData;

figure; plot(time_array, final.(animalID_to_use).(session_to_display).CNMFe_data.C_raw(neuron_to_display, :))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
% xline(BehavData.choiceTime(BehavData.shock == 1), '--y')
% add a plot for another cell, e.g. one that is not from the same epoch
% type
hold on; plot(time_array, final.(animalID_to_use).(session_to_display).CNMFe_data.C_raw(3, :))
