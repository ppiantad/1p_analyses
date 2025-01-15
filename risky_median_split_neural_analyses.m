% load 10x dataset or similar
% run generate_behav_figs_RDT.m to get risk_table. there is a column that
% contains a median split based on the % large reward choice across the 3
% blocks, use this to index into the variables created by
% block_wise_changes_v1.m

lost_all_risky = [lost_prechoice_sum(risk_table.risky == 1), lost_postchoice_sum(risk_table.risky == 1), lost_collection_sum(risk_table.risky == 1)]
conserved_all_risky = [conserved_prechoice_sum(risk_table.risky == 1), conserved_postchoice_sum(risk_table.risky == 1), conserved_collection_sum(risk_table.risky == 1)]
remapped_all_risky = [remapped_prechoice_sum(risk_table.risky == 1), remapped_postchoice_sum(risk_table.risky == 1), remapped_collection_sum(risk_table.risky == 1)]


lost_all_not_risky = [lost_prechoice_sum(risk_table.risky == 0), lost_postchoice_sum(risk_table.risky == 0), lost_collection_sum(risk_table.risky == 0)]
conserved_all_not_risky = [conserved_prechoice_sum(risk_table.risky == 0), conserved_postchoice_sum(risk_table.risky == 0), conserved_collection_sum(risk_table.risky == 0)]
remapped_all_not_risky = [remapped_prechoice_sum(risk_table.risky == 0), remapped_postchoice_sum(risk_table.risky == 0), remapped_collection_sum(risk_table.risky == 0)]




%% for plotting changes on a donut (specific) and pie (broad) charts
all_conserved_sum_risky = sum(conserved_all_risky)
all_lost_sum_risky = sum(lost_all_risky)
all_remapped_sum_risky = sum(remapped_all_risky)


total_neurons_risky = sum(neurons_per_mouse(risk_table.risky == 1))

remaining_neurons = total_neurons_risky - (all_conserved_sum_risky + all_lost_sum_risky +all_remapped_sum_risky);

figure;
piechart([all_conserved_sum_risky/total_neurons_risky, all_lost_sum_risky/total_neurons_risky, all_remapped_sum_risky/total_neurons_risky, remaining_neurons/total_neurons_risky])

figure;
donutchart([lost_prechoice_sum(risk_table.risky == 1), lost_all/neuron_num, remapped_all/neuron_num, remaining_neurons/neuron_num])

% all_neurons = [conserved_sum; lost_sum; remapped_sum]
%all_neurons_2 = [conserved_sum; lost_sum; remapped_sum]
%% for plotting changes on a donut (specific) and pie (broad) charts
all_conserved_sum_not_risky = sum(conserved_all_not_risky)
all_lost_sum_not_risky = sum(lost_all_not_risky)
all_remapped_sum_not_risky = sum(remapped_all_not_risky)

total_neurons_not_risky = sum(neurons_per_mouse(risk_table.risky == 0))

remaining_neurons = total_neurons_not_risky - (all_conserved_sum_not_risky + all_lost_sum_not_risky +all_remapped_sum_not_risky);

figure;
piechart([all_conserved_sum_not_risky/total_neurons_not_risky, all_lost_sum_not_risky/total_neurons_not_risky, all_remapped_sum_not_risky/total_neurons_not_risky, remaining_neurons/total_neurons_not_risky])

% figure;
% donutchart([conserved_all/neuron_num, lost_all/neuron_num, remapped_all/neuron_num, remaining_neurons/neuron_num])

% all_neurons = [conserved_sum; lost_sum; remapped_sum]
%all_neurons_2 = [conserved_sum; lost_sum; remapped_sum]


%% z test of proportions
% Variables
n1 = total_neurons_not_risky ; % Total sample size for "not risky"
n2 = total_neurons_risky ;     % Total sample size for "risky"

x1 = sum(remapped_all_not_risky); % Number of remapped neurons in "not risky"
x2 = sum(remapped_all_risky);     % Number of remapped neurons in "risky"

% Proportions
p1 = x1 / n1;
p2 = x2 / n2;

% Combined proportion
p_combined = (x1 + x2) / (n1 + n2);

% Z-statistic
z = (p1 - p2) / sqrt(p_combined * (1 - p_combined) * (1/n1 + 1/n2));

% Two-tailed p-value
p_value = 2 * (1 - normcdf(abs(z)));

% Display results
fprintf('Z-statistic: %.4f\n', z);
fprintf('P-value: %.4f\n', p_value);

if p_value < 0.05
    disp('The proportions are significantly different (p < 0.05).');
else
    disp('The proportions are not significantly different (p >= 0.05).');
end
