%% Decoding (with assistance from ChatGPT
start_time = 0; % sub-window start time
end_time = 3; % sub-window end time
ts1 = (-10:.1:10);

time_window = ts1 >= start_time & ts1 <= end_time;





% run acces_risk_inscopix_v2 on REW = 1.2
neuron_mean_large = neuron_mean;

%%
%run acces_risk_inscopix_v2 on REW = 0.3
neuron_mean_all = [neuron_mean_large; neuron_mean]
labels = [ones(size(neuron_mean_all,1)/2,1); zeros(size(neuron_mean_all, 1)/2,1)];

decoding_window = neuron_mean_all(:, time_window);

%%
cv = cvpartition(size(decoding_window, 1), 'HoldOut', 0.2);
idx_train = training(cv);
idx_test = test(cv);

%%
mdl = fitcsvm(decoding_window(idx_train,:), labels(idx_train));

%%
labels_pred = predict(mdl, decoding_window(idx_test,:));

%%
accuracy = sum(labels_pred == labels(idx_test)) / numel(labels_pred);