function tb = read_analysis_table(filename, table_name)
%READ_ANALYSIS_TABLE Read and validate an analysis-ready participant table.
if ~exist(filename, 'file')
    error('Missing %s input: %s', table_name, filename);
end
tb = readtable(filename, 'TreatAsMissing', {'NA','NaN',''});
validate_unique_ids(tb, table_name);
end
