% Load your DEG file
DEG_Table = readtable('PL-Nac_DEG (2).csv', 'PreserveVariableNames', true);

% Load the Ensembl-to-GeneName mapping file
Mapping_Table = readtable('MGI_Gene_Model_Coord.csv'); % Ensure this file has columns: 'EnsemblID' and 'GeneName'

% Create a dictionary for quick lookup
GeneMap = containers.Map(Mapping_Table.x11_EnsemblGeneId, Mapping_Table.x3_MarkerSymbol);

% Match gene names
numGenes = height(DEG_Table);
geneNames = cell(numGenes, 1);

for i = 1:numGenes
    id = DEG_Table{i, 1}{1}; % Extract Ensembl ID
    if isKey(GeneMap, id)
        geneNames{i} = GeneMap(id);
    else
        geneNames{i} = 'Unknown'; % If no match is found
    end
end

% Add the new column to the table
DEG_Table.GeneName = geneNames;

% Save the updated table
writetable(DEG_Table, 'PL-NAc_DEG_with_GeneNames.csv');

%%
DEG_Table.negLog10P = -log10(DEG_Table.padj);

log2FC = DEG_Table.log2FoldChange;
negLog10P = DEG_Table.negLog10P;
geneNames = DEG_Table.GeneName;

%%
% Define genes to annotate
genesToLabel = ["Lmo3", "Pcdh17", "Cd4", "Pnck"];

% Find indices of these genes
idx = ismember(geneNames, genesToLabel);

% Plot Volcano Plot
figure;
scatter(log2FC, negLog10P, 50, 'b', 'filled'); % Plot all points
hold on;
scatter(log2FC(idx), negLog10P(idx), 50, 'r', 'filled'); % Highlight selected genes

% Add annotations
text(log2FC(idx) + 0.1, negLog10P(idx), geneNames(idx), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');

% Labels and title
xlabel('Log_2 Fold Change');
ylabel('-log_{10}(P-value)');
title('Volcano Plot with Selected Gene Labels');
grid on;
legend({'Non-Significant', 'Highlighted Genes'}, 'Location', 'best');
hold off;