figure; imagesc(ts1, [], zall_mouse{2, 3}{1, 6})
figure; plot(ts1, mean(zall_mouse{2, 3}{1, 6}))


num_samples = size(final.BLA_Insc_25.Pre_RDT_RM.CNMFe_data.C_raw, 2)
% sampling_frequency = (final.(select_mouse).(first_session).choiceTime.uv.dt)*100;
sampling_frequency = (final.BLA_Insc_25.Pre_RDT_RM.uv.dt)*100;
% Create a time array
time_array = (0:(num_samples-1)) / sampling_frequency;

BehavData = final.BLA_Insc_25.Pre_RDT_RM.uv.BehavData;

figure; plot(time_array, final.BLA_Insc_25.Pre_RDT_RM.CNMFe_data.C_raw(6, :))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
% xline(BehavData.choiceTime(BehavData.shock == 1), '--y')
