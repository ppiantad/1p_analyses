for ii = 1:size(animalIDs,1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
        % only use rewarded trials for this, otherwise things get wonky
        [BehavData,trials,varargin]=TrialFilter_test(BehavData,'OMITALL', 0, 'BLANK_TOUCH', 0); 
        block_1 = [BehavData.stTime(BehavData.Block == 1) BehavData.collectionTime(BehavData.Block == 1)]; 
        block_1_mouse(ii,:) = [block_1(1, 1) block_1(end, 2)];
        block_2 = [BehavData.stTime(BehavData.Block == 2) BehavData.collectionTime(BehavData.Block == 2)]; 
        block_2_mouse(ii,:) = [block_2(1, 1) block_2(end, 2)];
        block_3 = [BehavData.stTime(BehavData.Block == 3) BehavData.collectionTime(BehavData.Block == 3)];
        block_3_mouse(ii,:) = [block_3(1, 1) block_3(end, 2)];

    end
end