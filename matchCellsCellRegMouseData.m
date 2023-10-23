data = mouseData;

%select your CellRegData folder containing all of the individual mouse
%folders.

curFolder = uigetdir;
cd(curFolder);
list = dir;
mask = [list.isdir];
list = list(mask);
list(1:2) = [];

clear mask

%get a list of each cellRegistered file in each folder, there should be 1
%per folder
for i = 1:length(list)
    cd(strcat(curFolder, '\',list(i).name))
    fileList = dir;
    
    
    idx = ~cellfun('isempty',strfind({fileList.name},'cellRegistered'));
    
    if sum(idx) > 1
        error('one of the folders has more than one cellRegistered mat file')
    elseif sum(idx) == 0
        error('one of the folders does not have a cellRegistered mat file')
    end
    
    x = find(idx, 1);
    
    cellFiles{i} = fileList(x).name;
    clear fileList idx x
end
clear i

%check if you have mice in the data struct for every registration file
if mod(length(data),length(cellFiles)) ~= 0
    error('The number of celRegistered files and mice in the data array are not matching properly')
end

%go through each registration entry and find the correct cells in the data
%array
for i = 1:length(list)
    
    cd(curFolder);
    cd(list(i).name);
    load(cellFiles{i});
    %This script assumes the first mouse is saline and the second is SKF
    %Find the matching entries in data
%     for t = 1:length(data)
% %         mouseID = data(t).mouseID
%         mouseID = data(t).mouseID
%         %         mouseIDs(t,:) = num2str(mouseID)
%     end
    mouseID = list(i).name;
    
    idx = ~cellfun('isempty',strfind({data.mouseID},list(i).name));
    %     idx = ~cellfun('isempty',strfind(mouseIDs,list(i).name));
    
    array = data(idx);
    
    treat = array(logical([array.day]));
    null = array(~(logical([array.day])));
    
    subArray(2) = null; %HAD TO CHANGE THIS FOR LEANDRA DATA. PREVIOUSLY NULL==1
    subArray(1) = treat;
    clear idx array treat null
    
    %align the datasets
    index = cell_registered_struct.cell_to_index_map;
    for k = 1:size(subArray,2)
        %         matchedArray(k).name = subArray(k).fileName;
        matchedArray(k).fileName = subArray(k).fileName;
        matchedArray(k).mouseID = subArray(k).mouseID;
        %         matchedArray(k).treatment = subArray(k).treatment;
        matchedArray(k).day = subArray(k).day;
        matchedArray(k).behav = subArray(k).behav;
        matchedArray(k).ca = subArray(k).ca;
        matchedArray(k).ca = rmfield(matchedArray(k).ca, 'caData');
        matchedArray(k).ca.caData = struct('rawTrace', {}, 'deconTrace', {}, 'eventTrace', {}, 'normalizedTraceZ', {}, 'normalizedTrace', {},'Coor',{} );
    end
    clear k
    
    %cycle through each cell and build the matched cell dataset
    counter = 1;
    for k = 1:size(index,1)
        if index(k,1) == 0
            continue
        elseif index(k,2) == 0
            continue
        end
        
        matchedArray(1).ca.caData(counter) = subArray(1).ca.caData(index(k,2)); %depending on order could be k,1 or k,2
        matchedArray(2).ca.caData(counter) = subArray(2).ca.caData(index(k,1));
        counter = counter +1;
    end
    clear counter k
    
    % for k = 1:size(subArray,2)
    %         matchedArray(k).numIcs = length(matchedArray(k).ca.data);
    % %         matchedArray(k).noldusBehavior = subArray(k).noldusBehavior;
    %         matchedArray(k).labJackBehavior = subArray(k).labJackBehavior;
    %         matchedArray(k).scopeBehavior = subArray(k).scopeBehavior;
    % %         matchedArray(k).offsets = subArray(k).offsets;
    %         matchedArray(k).discreteBehavior = subArray(k).discreteBehavior;
    % %         matchedArray(k).binaryTraces = subArray(k).binaryTraces;
    % end
    
    %add the data to the new mouse array
    if exist('mouseDataMatched') == 0
        mouseDataMatched = matchedArray;
    else
        mouseDataMatched = [mouseDataMatched, matchedArray];
    end
    clear matchedArray k counter index subArray cell_registered_struct
end

clear cellFiles curFolder data i list


