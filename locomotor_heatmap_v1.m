% run locomotor_trajectories_not_downsampled_data.m first! 

% Load data (replace these with your actual data)
% load('mouse_data.mat'); % assuming your data is stored in a .mat file
% X, Y coordinates sampled at 30 FPS
% time variable


block_1_length = BehavData.TrialPossible(BehavData.Block == 1);
block_1_ind = SLEAP_data.idx_time > block_1_length(1) & SLEAP_data.idx_time < block_1_length(end);

X_heatmap = SLEAP_data.x_pix(block_1_ind);
Y_heatmap = SLEAP_data.y_pix(block_1_ind);


% Create 2D histogram
nbins = 30; % number of bins for the histogram
histogram2D = hist3([X_heatmap Y_heatmap], [nbins nbins]);

% Normalize histogram
total_time = max(SLEAP_data.idx_time(block_1_ind)) - min(SLEAP_data.idx_time(block_1_ind)); % total time
histogram2D_normalized = histogram2D / total_time;

% Visualize heatmap
figure;
imagesc(histogram2D_normalized);
colormap('jet'); % Choose colormap (e.g., hot for brighter colors)
colorbar; % Add colorbar to show scale
title('Mouse Heatmap');
xlabel('X Coordinate');
ylabel('Y Coordinate');

%%





block_3_length = BehavData.TrialPossible(BehavData.Block == 3);
block_3_ind = SLEAP_data.idx_time > block_3_length(1) & SLEAP_data.idx_time < block_3_length(end);

X_heatmap = SLEAP_data.x_pix(block_3_ind);
Y_heatmap = SLEAP_data.y_pix(block_3_ind);


% Create 2D histogram
nbins = 30; % number of bins for the histogram
histogram2D = hist3([X_heatmap Y_heatmap], [nbins nbins]);

% Normalize histogram
total_time = max(SLEAP_data.idx_time(block_3_ind)) - min(SLEAP_data.idx_time(block_3_ind)); % total time
histogram2D_normalized = histogram2D / total_time;

% Visualize heatmap
figure;
imagesc(histogram2D_normalized);
colormap('jet'); % Choose colormap (e.g., hot for brighter colors)
colorbar; % Add colorbar to show scale
title('Mouse Heatmap');
xlabel('X Coordinate');
ylabel('Y Coordinate');

%%

block_3_length = BehavData.TrialPossible(BehavData.Block == 3);
block_3_ind = SLEAP_data.idx_time > block_3_length(1) & SLEAP_data.idx_time < block_3_length(end);

X_heatmap_block_3 = SLEAP_data.x_pix(block_3_ind);
Y_heatmap_block_3 = SLEAP_data.y_pix(block_3_ind);


xyPos = [X_heatmap_block_3,Y_heatmap_block_3]; %combine x and y data into a single two column matrix if they are not already

xyHist = hist3(xyPos,[100 100]); %generate bivariate histogram of xy data, gridsize 100x100

xyHistFilt = imgaussfilt(xyHist,4); %gaussian filter data, filtering value can be changed depending on data input
figure;
% Plot heatmap of filtered data
imagesc(xyHistFilt);
colormap hot


%%
block_1_length = BehavData.TrialPossible(BehavData.Block == 1);
block_1_ind = SLEAP_data.idx_time > block_1_length(1) & SLEAP_data.idx_time < block_1_length(end);

X_heatmap_block_1 = SLEAP_data.x_pix(block_1_ind);
Y_heatmap_block_1 = SLEAP_data.y_pix(block_1_ind);



block_3_length = BehavData.TrialPossible(BehavData.Block == 3);
block_3_ind = SLEAP_data.idx_time > block_3_length(1) & SLEAP_data.idx_time < block_3_length(end);

X_heatmap_block_3 = SLEAP_data.x_pix(block_3_ind);
Y_heatmap_block_3 = SLEAP_data.y_pix(block_3_ind);





% Calculate histograms for block 1 and block 3
xyPos_block_1 = [X_heatmap_block_1, Y_heatmap_block_1];
xyHist_block_1 = hist3(xyPos_block_1, [150 150]);

xyPos_block_3 = [X_heatmap_block_3, Y_heatmap_block_3];
xyHist_block_3 = hist3(xyPos_block_3, [150 150]);

% Gaussian filter the histograms
xyHistFilt_block_1 = imgaussfilt(xyHist_block_1, 4);
xyHistFilt_block_3 = imgaussfilt(xyHist_block_3, 4);

% Rotate the data by 90 degrees
xyHistFilt_block_1_rotated = rot90(xyHistFilt_block_1);
xyHistFilt_block_3_rotated = rot90(xyHistFilt_block_3);

% Determine the global min and max for both histograms
global_min = min(min(min(xyHistFilt_block_1_rotated)), min(min(xyHistFilt_block_3_rotated)));
global_max = max(max(max(xyHistFilt_block_1_rotated)), max(max(xyHistFilt_block_3_rotated)));

% Plot subplots with shared color scale
figure;

% Plot heatmap for block 1
subplot(1, 2, 1);
imagesc(xyHistFilt_block_1_rotated);
colormap hot;
caxis([global_min, global_max-18]); % Set color limits
title('Block 1 Heatmap (Rotated)');
colorbar;

% Plot heatmap for block 3
subplot(1, 2, 2);
imagesc(xyHistFilt_block_3_rotated);
colormap hot;
caxis([global_min, global_max-18]); % Set color limits
title('Block 3 Heatmap (Rotated)');
colorbar;

% Adjust subplot layout
sgtitle('Comparison of Block 1 and Block 3 Heatmaps (Rotated)');

%%

