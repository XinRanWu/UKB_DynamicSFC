function value = load_mat_variable(filename, variable_name)
%LOAD_MAT_VARIABLE Load one explicitly named variable from a MAT file.
if ~exist(filename, 'file')
    error('MAT file not found: %s', filename);
end
content = load(filename, variable_name);
if ~isfield(content, variable_name)
    error('%s does not contain variable %s.', filename, variable_name);
end
value = content.(variable_name);
end
