% Assuming you already have the structure 'final' loaded into the workspace

% Extract the substructure final.BLA_Insc_24
substruct = final.BLA_Insc_24.RDT_D1;

% Save the substructure and its nested contents to a .mat file
save('BLA_Insc_24_RDT_D1_Ca_data.mat', '-struct', 'substruct');

%%

% Extract the substructure final.BLA_Insc_24
substruct = final_SLEAP.BLA_Insc_24.RDT_D1;

% Save the substructure and its nested contents to a .mat file
save('BLA_Insc_24_RDT_D1_SLEAP_data.mat', '-struct', 'substruct');
