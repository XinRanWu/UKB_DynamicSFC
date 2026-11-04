function out = innerjoin_many(varargin)
%INNERJOIN_MANY Inner-join one-row-per-eid tables with duplicate checks.
if nargin < 1
    error('At least one table is required.');
end
out = varargin{1};
validate_unique_ids(out, 'input table 1');
for i = 2:nargin
    validate_unique_ids(varargin{i}, sprintf('input table %d', i));
    overlapping = intersect(out.Properties.VariableNames, ...
        varargin{i}.Properties.VariableNames, 'stable');
    overlapping(strcmp(overlapping, 'eid')) = [];
    if ~isempty(overlapping)
        error('Duplicate non-eid variables before join: %s', ...
            strjoin(overlapping, ', '));
    end
    out = innerjoin(out, varargin{i}, 'Keys', 'eid');
end
end
