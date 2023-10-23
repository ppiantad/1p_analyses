%% behavioral raster script. this script will generate pseudo licking data
% to apply to your own behavioral data, you will need a binarized trial x behavior matrix
% where 0s = no behavior and 1 = behavior (e.g. lick) occurred
% the script generates fake trial aligned behavioral data
%% generate fake lick data
preLicks = (rand(50,200)>0.95); %generate sparse licking pre event
eventLicks = (rand(50,200)>0.8); %generate dense licking during event
postLicks = (rand(50,200)>0.95); %generate sparse licking post event
%% combine all licks
allLicks = horzcat(preLicks,eventLicks,postLicks); %concat pseudolicks
allLicks = double(allLicks); %convert from logical to an array
%% convert 0s to NaNs for plotting, if you dont do this zeros will be plotted. start here if you have trial x behavior matrix of 0s and 1s ready
allLicks(allLicks==0) = NaN; %allLicks would be your binarized behavioral events (trial x behavior)
%% plot data, modify this if you have your behavioral matrix read to go
figure
hold
count = 0;
for x = 1:size(allLicks)
    evtWin = -19.9:0.10:40; %you will need to change this to fit your data
    plot(evtWin,allLicks(x,:) + count,'|k') % | symbol for plotting
    count = count + 1;
end
xline(0,'-'); %mark where event occurred
set(gca,'FontName','Arial','FontSize',16)
set(gcf, 'Position',  [100, 100, 300, 600])
xlabel('Time from event (s)','FontSize',22)
ylabel('Trial number','FontSize',22)