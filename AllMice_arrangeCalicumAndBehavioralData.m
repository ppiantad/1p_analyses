
% clear all; close all; clc;
% %%
% cd('C:\Users\1\Desktop\temp files\ca Event Analysis')
% load('mouseArrayData')

data = mouseData;
%% select mouse to run analyses on

for m = 1:length(data)

%% grab selected data (user modifiable code)
tempBehTS =  mouseData(m).discreteBehavior.Grooming.timeStamp;          %get behavior timestamps for selected event 
tempBehStart = mouseData(m).discreteBehavior.Grooming.startIdx;         %direct transfer of logical index for starts of behavior
tempBehStop = mouseData(m).discreteBehavior.Grooming.stopIdx;           %direct transfer of logical index for stops of behavior
%
behav.GroomingStart = tempBehTS(tempBehStart);                          %use logical index to select the starts of behaviordirect transfer of 
behav.GroomingStop = tempBehTS(tempBehStop);                            %use logical index to select the stop times for a behavior
%
clear temp*


% %%
% tempBehTS =  mouseData(m).discreteBehavior.FaceGrooming.timeStamp;
% tempBehStart = mouseData(m).discreteBehavior.FaceGrooming.startIdx;
% tempBehStop = mouseData(m).discreteBehavior.FaceGrooming.stopIdx;
% %
% behav.faceGroomingStart = tempBehTS(tempBehStart);
% behav.faceGroomingStop = tempBehTS(tempBehStop);
% %
% clear temp*



%%
ca.time = mouseData(m).caData.eventTraceInfo.timeVector;                    %direct transfer of calcium trace's time vector
%%
for u = 1:size(mouseData(m).caData.data,2)                                  %for each unit
    ca.events(u,:) = mouseData(m).caData.data(u).eventTrace;                %direct transfer of calcium event indicator (not binary; includes magnitude)
    idx = ca.events(u,:) ~= 0;                                              %logical index events (non-zero time bins)
    ca.eventTS{u} = ca.time(idx);                                           %time of calcium events (indexed time values)
    %%
    ca.traces(u,:) = mouseData(m).caData.data(u).rawCaTrace;                   %direct transfer of calcium trace
end
% clear u idx mouseData m

clear u idx

processed(m).name = data(m).name;
processed(m).mouseID = data(m).mouseID;
processed(m).treatment = data(m).treatment;
processed(m).numIcs = data(m).numIcs;
processed(m).behav = behav;
processed(m).ca = ca;
clear behav ca
end

clear m data
