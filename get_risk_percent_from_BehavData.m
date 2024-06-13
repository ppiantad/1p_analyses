% iter = iter+1;
% neuron_num = 0;
risk_table = table;
for ii = 1:size(animalIDs,1) % 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
        % only use rewarded trials for this, otherwise things get wonky
        [BehavData,trials,varargin]=TrialFilter(BehavData,'ALL', 1); 
        block_1_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 1 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 1); 
        % block_1_mouse(ii,:) = [block_1(1, 1) block_1(end, 2)];
        block_2_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 2 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 2);  
        % block_2_mouse(ii,:) = [block_2(1, 1) block_2(end, 2)];
        block_3_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 3 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3);  
        % block_3_mouse(ii,:) = [block_3(1, 1) block_3(end, 2)];
        lose_shift_percent = sum(BehavData.lose_shift == 1)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3);
        lose_omit_percent = sum(BehavData.lose_omit == 1)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3);
        lose_stay_percent = sum(BehavData.lose_stay == 1)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3);
        large_aborts = []; %[] sum(BehavData.type_binary == 1)
    end
    risk_table(ii,:) = array2table([block_1_large_choice_percent, block_2_large_choice_percent, block_3_large_choice_percent, large_aborts, lose_shift_percent, lose_omit_percent, lose_stay_percent]);
end

% some mice have NaNs if they didn't make it to this trial block. replace
% the NaNs with 0 because not making it to the trial block is basically
% being as risk averse as possible. 
risk_table{:, :}(isnan(risk_table{:, :})) = 0;
row_means = nanmean(risk_table{:, 1:3}, 2);
risk_table.Mean_1_to_3 = row_means;
riskiness = risk_table.Mean_1_to_3;
aborts = risk_table.Var4;

