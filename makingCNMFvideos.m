%%Load in preprocessed data file then run lines indivually below

%%
% t_begin = 3000; %start frame you want in video
% t_end = 4000; %end frame you want in video

tmp_range = [1,1000];

%tmp_range = [t_begin, min(t_begin+100*kt-1, t_end)];
%tmp_range = [t_begin, min(t_begin+100*kt-1, t_end)];

Y = neuron.load_patch_data([], tmp_range);
Ybg = neuron.reconstruct_background(tmp_range);
Ysignal = double(Y) - Ybg; %this works!
%denoised data
%%
% Ydenoised = neuron.reshape(neuron.A*neuron.C,2);

% neuron.playMovie(Y,[1000 2500],'gray','Raw'); %adding 'Raw' will save the
% %video with that name (can be called whatever you want)

% neuron.playMovie(Y,[300 500],'gray'); %raw signal
% 
% neuron.playMovie(Ybg,[200 1000],'gray'); %background
% 
% % test = normalize(Ysignal); %normalize Ysignal
% 
neuron.playMovie(Ysignal,[-10 40],'parula','bckground_subtract_heatmap');

% neuron.playMovie(Ydenoised,[],'parula');
% 
% test = correlation_image(Ysignal); %correlation image of signal
% imagesc(test)
% colormap('parula')
% 
% Coor = neuron.show_contours(0.6);
% 
% % %% for tone info
% % range_ac = [-5, 70]
% % 
% % for tt=t_begin:1:t_end
% %     m = tt;
% % imagesc(double(Y(:, :, m))-Ybg(:, :, m), range_ac); hold on;
% % %Ysignal = double(Y(:, :, m))-Ybg(:, :, m);
% % text(1, 10, sprintf('Time = %.2f', tt/10), 'fontsize', 15, 'color', 'w');
% % %text(1, 10, sprintf('Time = %.2f', t/obj.Fs), 'fontsize', 15, 'color', 'w');
% % drawnow();
% % colormap('gray')
% % end


