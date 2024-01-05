%% Run eventRelatedActivity first for whatever events you want to identify responsive neurons for

select_mouse = 'BLA_Insc_25';

select_mouse_index = find(strcmp(animalIDs, select_mouse));

figure; plot(ts1, mean(neuron_mean_mouse{select_mouse_index, 2}(respClass_mouse.(select_mouse).Pre_RDT_RM.choiceTime.Outcome_Minus_4to0.REW_Large.activated== 1, :)))
hold on; plot(ts1, mean(neuron_mean_mouse{select_mouse_index, 3}(respClass_mouse.(select_mouse).Pre_RDT_RM.collectionTime.Outcome_1to3.REW_Large.activated == 1, :)))


%% organize data - lots of this is manual due to need to specify which events have been previously filtered


caTraceTrials = final.(select_mouse).Pre_RDT_RM.CNMFe_data.C_raw(respClass_mouse.(select_mouse).Pre_RDT_RM.choiceTime.Outcome_Minus_4to0.REW_Large.activated== 1, :);

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

caTraceTrials = final.(select_mouse).Pre_RDT_RM.CNMFe_data.C_raw(respClass_mouse.(select_mouse).Pre_RDT_RM.collectionTime.Outcome_1to3.REW_Large.activated== 1, :);

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
%%
% in order to trim off excess (because calcium recording starts before
% behavior), you'll need the start time. unfortunately I don't save that
% variable in the "final" struct, but you can get it from the adjustment to
% some of the columns of BehavData. e.g., see below - but make sure to
% update the session etc as necessary! 

BehavData = final.(select_mouse).Pre_RDT_RM.choiceTime.uv.BehavData;
% because the first trial possible is ALWAYS 60 seconds after ABET is
% issued, you can determine what adjustment has been made to this column
% (adding time to account for calcium recording starting first) by
% subtracting off 60 from the first element
stTime = BehavData.TrialPossible(1)-60; 

%%
num_samples = size(zall_second_event, 2)
sampling_frequency = (final.(select_mouse).Pre_RDT_RM.choiceTime.uv.dt)*100;
% Create a time array
time_array = (0:(num_samples-1)) / sampling_frequency;

%probably will want to filter these data to smooth them? let's use the
%inferred spike data for now
%chosing a random mouse
figure; plot(time_array, mean(zall_first_event))
hold on; plot(time_array, mean(zall_second_event))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2), '--k')




