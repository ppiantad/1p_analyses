% Helper function to intersect multiple tables
function combined_table = intersect_tables(tables)
    if isempty(tables)
        combined_table = [];
        return;
    end
    
    combined_table = tables{1};
    for i = 2:length(tables)
        combined_table = intersect(combined_table, tables{i}, 'rows');
    end
end