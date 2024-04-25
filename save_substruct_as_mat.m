% Assuming you already have the structure 'final' loaded into the workspace

% Extract the substructure final.BLA_Insc_24
substruct = final.BLA_Insc_24;

% Create a new structure with the same hierarchy
new_struct = struct('BLA_Insc_24', substruct);

% Save the new structure to a .mat file
save('substruct_file.mat', '-struct', 'new_struct');

%%

% Extract the substructure final.BLA_Insc_24
substruct = final_SLEAP.BLA_Insc_24.RDT_D1;

% Save the substructure and its nested contents to a .mat file
save('BLA_Insc_24_SLEAP_data.mat', '-struct', 'substruct');