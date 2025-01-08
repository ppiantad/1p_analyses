%WRITTEN BY BRITT CHAMBERLAIN 11232024 

%EDITED 01062025 

%We now have the data projected into PC space fitted to a reference session
%here, each matrix is 15x41: 15 = the top 15 principal components, 41 = 41
%timepoints (from -2s to 2s, sampled at 10Hz) 

%we want to calculate the Euclidean distance at every point in time: the
%equation for euclidean distance (over the first 3 PCs) is: 

    %d(p,q) = sqrt((q1 - p1)^2 + (q2 - p2)^2 + q3 - p3)^2)

    clc;
clear;
%OPEN FILE GENERATED IN PYTHON: 
[meta_import1, meta_folder1] = uigetfile('*.*','All Files (*.*)', 'MultiSelect','on');

meta_import1 = meta_import1'; 
name = convertCharsToStrings(meta_import1); 
load(meta_import1); 

%example: compare S1 and S7: 
s7 = s7pcaData(1:3,:); 
s1 = projS1(1:3,:); 

eucD1to7= [];
%here, p = s1 and q = s7
for i = 1:size(s1,2)
    p1 = s1(1,i);
    p2 = s1(2,i);
    p3 = s1(3,i);

    q1 = s7(1,i);
    q2 = s7(2,i);
    q3 = s7(3,i); 

    eucDtemp = sqrt((q1-p1).^2 + (q2-p2).^2 + (q3-p3).^2);
eucD1to7(1,i) = eucDtemp;
end 


%compare this to S6 vs S7: 
s6 = projS6(1:3,:);

eucD6to7 = [];
%now, p = s6 and q = s7
for i = 1:size(s6,2) 
    p1 = s6(1,i);
    p2 = s6(2,i);
    p3 = s6(3,i); 

    q1 = s7(1,i);
    q2 = s7(2,i);
    q3 = s7(3,i); 

    eucDtemp = sqrt((q1-p1).^2 + (q2-p2).^2 + (q3-p3).^2)

    eucD6to7(1,i) = eucDtemp;
end 

%I'm mostly interested in the distance between the trajectories from 500ms
%before avoid response through 1sec after: I can compare the euclidean
%distances at these points in time / get the average euclidean distance
%during this period of time: 
restrictWin1 = eucD1to7(1,15:31);
restrictWin6= eucD6to7(1,15:31); 
avgEucDist_1to7 = nanmean(restrictWin1);
avgEucDist_6to7 = nanmean(restrictWin6); 
