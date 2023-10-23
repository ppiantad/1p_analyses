load('BLA_Risk_Data_struct.mat')


session_to_analyze = 'RDT_D1';
epoc_to_align = 'stTime';
event_to_analyze = {'REW',1.2};

ts1 = (0:.1:20);
%%

animalIDs = (fieldnames(final));

for ii = 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
   
    if isfield(final.(currentanimal), session_to_analyze)
        [data,trials,~] = TrialFilter(final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.BehavData,'REW',1.2);
        trials = cell2mat(trials);
        for qq = 1:size(final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials,2)
            cell_zscores{qq,:} = final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(qq).zall(trials,(1:length(ts1)));
            cell_zscore_mean(qq,:) = mean(final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(qq).zall(trials,(1:length(ts1))));


            cell_raw{qq,:} = final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(qq).caTraces(trials,(1:length(ts1)));
            cell_raw_mean(qq,:) = mean(final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(qq).caTraces(trials,(1:length(ts1))));


%             kk = 1;
%             for kk = 1:size(final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(kk).zall,1)
%                 zall_cell{ii,:} =  
%             end
        end
    elseif ~isfield(final.(currentanimal), session_to_analyze)
        

    end
end