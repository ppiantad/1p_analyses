%% Run kmeans on all mice together
[coeff,score,latent,tsquared,explained] = pca(neuron_mean); %PCA to see how many clusters

%plot(explained(1:10,:)) %scree plot

[idx,C,sumdist3] = kmeans(neuron_mean,4,'Distance','correlation','Display','final', 'Replicates', 200,'Start','uniform');