%WRITTEN BY BRITT CHAMBERLAIN 11242024 

%population-level stabilization of avoidance ensembles: 
%here, we want to look at the activity of matched cells during an event of
%interest (in my case, avoidance responding) across multiple sessions and
%assess how correlated population activity is from day to day. 

%here, I am looking at neural activity from the very first session (the first few avoid responses a mouse is making) 
%all the way through very well learned behavior ("optimal performance") on the final (7th) session 

%I'd predict that population representation of avoidance would be most similar between late sessions
%as compared to the first session vs. last session. 

%I'm assessing this by taking data from one session and performing PCA to
%reduce dimensionality
    %then, I will PROJECT population data from other sessions into this
    %EXISTING PC space

    %at every point in time, we can assess the Euclidean distance between
    %our reference session and each individual session. Smaller Euclidean
    %distance = more similar. 

%Here, I have a file with the neural activity and time bins of behavioral
%events from all matched cells ("exampleFile_allSessions_matchedData.mat")

%FOR EACH SESSION, EXTRACT: 
    % 1) the timestamp of all avoid responses 
    % 2) the session normalized trace for all matched cells 

    %use these to get the trial-averaged response for each neuron 
    %save this in a structure to then perform PCA 

clc;
clear;
%SELECT FILE OF INTEREST 
[meta_import1, meta_folder1] = uigetfile('*.*','All Files (*.*)', 'MultiSelect','on');

meta_import1 = meta_import1'; 
name = convertCharsToStrings(meta_import1); 
load(meta_import1); 
%% 
fn = fieldnames(matchData); 
trajectoryInputs = struct; 

for t = 1:numel(fn) 

    fn1 = matchData.(fn{t}); 
    avoidIdx = fn1.eventStamps.avoidResponse; %extract event indices
    numTrials = length(avoidIdx); 

    normCells = fn1.cellTrace.normTrace; %extract session normalized trace for each cell
    numCells = size(normCells,2); 

    %determine the window of time around your event that you want to look
    %at (ex: 2s prior to 2s after) 
    evPre = 20; %2 seconds prior to event (avoid response) occurring
    evPost = 20; %2seconds after event occurring


    trialAvgResponse = [];
    for i = 1:numCells 
        tempEv = [];
        for j = 1:numTrials
            evWin = normCells((avoidIdx(j)-evPre):(avoidIdx(j)+evPost),i); 
            tempEv = [tempEv,evWin];
        end 
        trialAvg_eachCell = mean(tempEv,2);
        trialAvgResponse = [trialAvgResponse,trialAvg_eachCell];
    end 
    trialAvgResponse = trialAvgResponse'; 

    %save this into the larger structure: 
    sessionNum = strcat('s',num2str(t));
    trajectoryInputs.(sessionNum) = trialAvgResponse;
end 

clearvars -except matchData fn trajectoryInputs name


cd avoidTrajectoryInputs; 

reducedName = extractBefore(name,'matchedData'); 
save(strcat(reducedName,'avoidTrajectoryInputs.mat'),'trajectoryInputs'); 

cd .. 
clc;
clear;

%%FROM HERE, WE SWITCH TO PYTHON TO PERFORM PCA 
    %WHY? I've found that python has better functions for
    %transforming/projecting data into PC space. It also "centers" the data
    %