for ii = 1:size(animalIDs,1) % 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
        [BehavData,trials,varargin]=TrialFilter_test(BehavData, 'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 3);

        y_outcomes = BehavData.bigSmall;
        x_latencies = BehavData.choTime2;


        % Fit a linear model
        mdl = fitlm(x_latencies, y_outcomes);

        % Display the model summary
        disp(mdl)

        % Extract the coefficient for latency
        latency_coef(ii) = mdl.Coefficients.Estimate(2);
        p_value(ii) = mdl.Coefficients.pValue(2);

        % % Display the results
        % fprintf('Coefficient for latency: %.4f\n', latency_coef);
        % fprintf('P-value for latency: %.4f\n', p_value);
    end
end
%%

for ii = 1:size(animalIDs,1) % 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
        [BehavData,trials,varargin]=TrialFilter_test(BehavData, 'OMITALL', 0, 'BLANK_TOUCH', 0);
        y_outcomes = BehavData.bigSmall;
        x_latencies = BehavData.choTime2;
        % Assuming y_outcomes is your vector of choices
        y_outcomes(y_outcomes == 0.3) = 0;
        y_outcomes(y_outcomes == 1.2) = 1;
        % Fit a generalized linear model (assuming binary choice 0 or 1)
        [b, dev, stats] = glmfit(x_latencies, y_outcomes, 'binomial', 'link', 'logit');

        % Extract the coefficient for latency
        latency_coef(ii) = b(2);
        p_value(ii) = stats.p(2);
        % 
        % % Display the results
        % fprintf('Coefficient for latency: %.4f\n', latency_coef);
        % fprintf('P-value for latency: %.4f\n', p_value);
    end
end