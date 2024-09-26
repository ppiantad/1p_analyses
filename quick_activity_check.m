animalIDs = (fieldnames(final));
select_mouse = 'BLA_Insc_40';

select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'RDT_D1';

BehavData = final.(select_mouse).(first_session).uv.BehavData;
% because the first trial possible is ALWAYS 60 seconds after ABET is
% issued, you can determine what adjustment has been made to this column
% (adding time to account for calcium recording starting first) by
% subtracting off 60 from the first element
stTime = BehavData.TrialPossible(1)-60; 

caTraceTrials = final.(select_mouse).(first_session).CNMFe_data.C_raw;

ca_data = final.(select_mouse).(first_session).CNMFe_data.C_raw;
ca_data_denoised = final.(select_mouse).(first_session).CNMFe_data.C;
ca_data_spike_inf = final.(select_mouse).(first_session).CNMFe_data.spike_prob;
ca_data_cnmfe_spike = full(final.(select_mouse).(first_session).CNMFe_data.S);
ca_data_norm = normalize(ca_data, 2);
ca_data_denoised_norm = normalize(ca_data_denoised, 2);
ca_data_spike_inf_norm = normalize(ca_data_spike_inf, 2);


num_samples = size(ca_data, 2)
% sampling_frequency = (final.(select_mouse).(first_session).choiceTime.uv.dt)*100;
sampling_frequency = (final.(select_mouse).(first_session).uv.dt)*100;

time_array = (0:(num_samples-1)) / sampling_frequency;

figure;
plot(time_array, ca_data_norm(5, :))
% hold on; plot(time_array, ca_data_denoised(1, :))
hold on; plot(time_array, ca_data_cnmfe_spike(5, :))
% hold on; plot(time_array, choiceVar)
% hold on; plot(time_array, ca_data_spike_inf_norm(13, :))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
xline(BehavData.choiceTime(BehavData.shock == 1), '--y')