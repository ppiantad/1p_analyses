iter = iter+1;
neuron_num = 0;
risk_table = table;
for ii = 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
        % only use rewarded trials for this, otherwise things get wonky
        [BehavData,trials,varargin]=TrialFilter(BehavData,'OMITALL', 0, 'BLANK_TOUCH', 0); 
        block_1_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 1 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 1); 
        % block_1_mouse(ii,:) = [block_1(1, 1) block_1(end, 2)];
        block_2_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 2 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 2);  
        % block_2_mouse(ii,:) = [block_2(1, 1) block_2(end, 2)];
        block_3_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 3 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3);  
        % block_3_mouse(ii,:) = [block_3(1, 1) block_3(end, 2)];


    end
    risk_table(ii,:) = array2table([block_1_large_choice_percent, block_2_large_choice_percent, block_3_large_choice_percent]);
end

risk_table(:, end+1) = nanmean(risk_table, 2);
riskiness = table2array(nanmean(risk_table, 2));