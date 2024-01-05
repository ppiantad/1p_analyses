% Assuming you have your original struct named "cellreg_struct"
% and the field you want to extract is "Pre_RDT_RM_vs_RDT_D1"

% Create the "final" structure
final = struct();

% Extract the "Pre_RDT_RM_vs_RDT_D1" field from cellreg_struct
preRDTField = cellreg_struct.RM_D1_vs_Pre_RDT_RM;
% preRDTField = cellreg_struct.Pre_RDT_RM_vs_RDT_D1;
% preRDTField = cellreg_struct.RDT_D1_vs_RDT_D2;
% preRDTField = cellreg_struct.RDT_D2_vs_RDT_D3;
% preRDTField = cellreg_struct.RDT_D1_vs_SHOCK_TEST;

% Copy the data from preRDTField and its subfields into final
final = preRDTField;

% If you want to copy all fields inside preRDTField to final, you can use a loop
% for example, if preRDTField is a struct itself:
% fields = fieldnames(preRDTField);
% for i = 1:length(fields)
%     final.(fields{i}) = preRDTField.(fields{i});
% end