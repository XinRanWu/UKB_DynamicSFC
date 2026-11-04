function validate_unique_ids(tb, table_name)
%VALIDATE_UNIQUE_IDS Validate a one-row-per-participant table.
if ~istable(tb) || ~ismember('eid', tb.Properties.VariableNames)
    error('%s must be a table containing eid.', table_name);
end
if any(ismissing(tb.eid))
    error('%s contains missing eid values.', table_name);
end
if numel(unique(tb.eid)) ~= height(tb)
    error('%s contains duplicated eid values.', table_name);
end
end
