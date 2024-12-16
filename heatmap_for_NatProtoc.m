%
load('heatmap_for_NatProtoc.mat')

%%
custom_colormap = [
    1, 1, 1; % white
    1, 0.9, 0.9;
    1, 0.8, 0.8;
    1, 0.7, 0.7;
    1, 0.6, 0.6;
    1, 0.5, 0.5;
    1, 0.4, 0.4;
    1, 0.3, 0.3;
    1, 0.2, 0.2;
    1, 0.1, 0.1;
    1, 0, 0;   % red
];


n = 256; 
custom_colormap = interp1(linspace(0, 1, size(custom_colormap, 1)), custom_colormap, linspace(0, 1, n));


figure('Position', [100, 100, 350, 600]); % [left, bottom, width, height]
hold on


imagesc(ts1, 1, neuron_mean_sorted);


colormap(custom_colormap);

clim([-1 1]);


c = colorbar('eastoutside');
set(c, 'YTick', clim); % 
ylim([1, neuron_num]);
xlim([-8 8]);

set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, neuron_num]);
% xline(0)

fontname('arial')
fontsize(20, 'points')
hold off;