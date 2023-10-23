function [animalNames, blockpaths, behavFiles, session, implant_side, large_rew_side, whichStreams, whichTTL] = read_mouse_data(root, which_region, which_sessions);

if which_region == '1'
    mouse_ID={
        'BLA-Insc-1';...
        'BLA-Insc-2';...
        'BLA-Insc-3';...
        'BLA-Insc-4';...
    
        }
    for i=1:size(mouse_ID)
        
RDT_D1_files={
    'BLA-INSC-6 05182021.csv','BLA-Insc-6_2021-05-18_RDT_D1_GPIO.csv','BLA-Insc-6_2021-05-18_RDT_D1.CNMF_final.mat';...
    'BLA-INSC-3','EAYEAYAEY','EYAYAE';...
         };

blockpaths={
    strcat(root,'\MATLAB\TDTbin2mat\Photometry\102 & 103\102_and_103-200113-131356')...
    strcat(root,'\MATLAB\TDTbin2mat\Photometry\102 & 103\102_and_103-200114-130456')...
    };

behavFiles={
    'RRD102 01132020.csv'...
    'RRD102 01142020.csv'...
    };

session = [1, 2];
% %Side of photometry implant, 1 = Left; 2 = Right
implant_side=[1, 1];
% %Which screen is the large rew side, 1 = Left; 2 = Right
large_rew_side=[1, 1];
% Sensor D = 12; Sensor C = 34
whichStreams=[34; 34];
% Old recordings have 2 TTLs, newer ones include "Tick" or "CAM1" which
% requires > 2. MUST CHECK data.streams.epocs for each mouse to make sure
% the correct TTL is selected for that particular session (the order is
% dependent on the # of mice recorded from on that day, etc)
whichTTL=[4; 4];

end
end
