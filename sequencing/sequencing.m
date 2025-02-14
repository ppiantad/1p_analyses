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
% genesToLabel = ["Lmo3", "Pcdh17", "Cd4", "Pnck"];

genesToLabel = ["Lmo3", "Pcdh17", "Cd4", "Car1", "Cenpm", "Ddc", "Slc18a3"];

% Find indices of these genes
idx = ismember(geneNames, genesToLabel);

% Define thresholds
log2FC_threshold = 0.5;
negLog10P_threshold = 1;

% Define categories based on thresholds
isRed = (negLog10P > negLog10P_threshold) & (log2FC > log2FC_threshold | log2FC < -log2FC_threshold); % Significant
isGray = (log2FC >= -log2FC_threshold & log2FC <= log2FC_threshold) & (negLog10P < negLog10P_threshold); % Low FC, Not Significant
isGreen = (log2FC > log2FC_threshold | log2FC < -log2FC_threshold) & (negLog10P < negLog10P_threshold); % High FC, Not Significant
isBlue = ~(isRed | isGray | isGreen); % Default color for other points
isBlack = idx; % Genes that should be labeled separately

% Create Volcano Plot
figure; hold on;

% Scatter plots with different colors
scatter(log2FC(isGray), negLog10P(isGray), 50, [0.5 0.5 0.5], 'filled'); % Grey points
scatter(log2FC(isGreen), negLog10P(isGreen), 50, 'g', 'filled'); % Green points
scatter(log2FC(isRed), negLog10P(isRed), 50, 'r', 'filled'); % Red points
scatter(log2FC(isBlue), negLog10P(isBlue), 50, 'b', 'filled'); % Default blue points

% Scatter plot for genes in genesToLabel (Black)
scatter(log2FC(isBlack), negLog10P(isBlack), 80, 'k', 'filled', 'd'); % Black diamonds for visibility

% Add vertical and horizontal threshold lines
xline(log2FC_threshold, 'k--', 'LineWidth', 1.5);  % +0.5 Log2FC
xline(-log2FC_threshold, 'k--', 'LineWidth', 1.5); % -0.5 Log2FC
yline(negLog10P_threshold, 'k--', 'LineWidth', 1.5); % -Log10P = 1.5

% Label selected genes
text(log2FC(isBlack) + 0.1, negLog10P(isBlack), geneNames(isBlack), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');

% Labels, Title, and Legend
xlabel('Log_2 Fold Change');
ylabel('-log_{10}(P-value)');
legend({'Low FC & Not Significant', 'High FC & Not Significant', 'Significant', 'Other', 'Genes of Interest'}, ...
    'Location', 'best');
hold off;
