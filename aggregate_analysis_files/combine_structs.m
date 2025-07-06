% % Load the two .mat files
% data1 = load('BLA-NAcShell_PdCO_behavior_data_04072025.mat');  % contains final_behavior
% data2 = load('BLA-NAcSh_behavior_only_01082025.mat');           % contains final_behavior
% 
% % Create a meta struct
% final_behavior = struct();
% 
% % Assign the loaded structs as subfields for clarity
% final_behavior.PdCO = data1.final_behavior;
% final_behavior.NAcSh = data2.final_behavior;
% 
% % Save the meta struct to a new .mat file
% save('meta_final_behavior.mat', 'final_behavior');

%%
% Load both .mat files
data1 = load('BLA_panneuronal_behavior_12302024.mat');  % final_behavior
data2 = load('females_behavior.mat');          % final_behavior

% Extract structs
fb1 = data1.final_behavior;
fb2 = data2.final_behavior;

% Combine fields into one struct
final_behavior = fb1;

% Add all fields from fb2 into final_behavior
fields_fb2 = fieldnames(fb2);
for i = 1:numel(fields_fb2)
    field = fields_fb2{i};
    final_behavior.(field) = fb2.(field);
end




%Save the merged struct
save('females_with_inscopix_males_behavior.mat', 'final_behavior');

%% for SLEAP data
% Load both .mat files
data1 = load('females_SLEAP.mat');  % final_SLEAP
data2 = load('BLA_panneuronal_SLEAP_data_2024_11_04_plus_XY_correction.mat');          % final_SLEAP

% Extract structs
fb1 = data1.final_SLEAP;
fb2 = data2.final_SLEAP;

% Combine fields into one struct
final_SLEAP = fb1;

% Add all fields from fb2 into final_behavior
fields_fb2 = fieldnames(fb2);
for i = 1:numel(fields_fb2)
    field = fields_fb2{i};
    final_SLEAP.(field) = fb2.(field);
end




%Save the merged struct
save('females_with_inscopix_males_SLEAP.mat', 'final_SLEAP');

