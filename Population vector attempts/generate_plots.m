figure; plot(time_array, prechoice_similarityOverTime(1, 1:size(time_array, 1)))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
if strcmp('shock',BehavData.Properties.VariableNames)
    xline(BehavData.choiceTime(BehavData.shock == 1), '--y')
end

%%
figure; plot(time_array, postchoice_similarityOverTime(1, 1:size(time_array, 1)))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
if strcmp('shock',BehavData.Properties.VariableNames)
    xline(BehavData.choiceTime(BehavData.shock == 1), '--y')
end

%%
figure; plot(time_array, consumption_similarityOverTime_array{1, 1}  (1, 1:size(time_array, 1)))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
if strcmp('shock',BehavData.Properties.VariableNames)
    xline(BehavData.choiceTime(BehavData.shock == 1), '--y')
end
hold on; plot(time_array, consumption_similarityOverTime_array{2, 1}  (1, 1:size(time_array, 1)))

%%

consumption_similarityOverTime_zscore = zscore(consumption_similarityOverTime, [], 2);

figure; plot(time_array, consumption_similarityOverTime_zscore(1, 1:size(time_array, 1)))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
if strcmp('shock',BehavData.Properties.VariableNames)
    xline(BehavData.choiceTime(BehavData.shock == 1), '--y')
end