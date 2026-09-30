function results = LMM_withFormula_Interact(tb, formula, x1_list, x2_list, y_list)
% Fit linear mixed models for predictor-by-moderator interactions.
% Main analysis files are scripts; this utility is a function because the same
% model engine is reused by demographic and genetic analyses.

if numel(x2_list) == 1
    x2_list = repmat(x2_list, size(x1_list));
elseif numel(x2_list) ~= numel(x1_list)
    error('x2_list must have one item or match x1_list in length.');
end

covariates = extract_formula_variables(formula);
n_models = numel(y_list) * numel(x1_list);
result_cells = cell(n_models, 1);

parfor model_index = 1:n_models
    [y_index, x_index] = ind2sub( ...
        [numel(y_list), numel(x1_list)], model_index);
    y_name = y_list{y_index};
    x1_name = x1_list{x_index};
    x2_name = x2_list{x_index};
    needed = unique([{y_name, x1_name, x2_name}, covariates], 'stable');
    model_table = tb(:, needed);

    numeric_names = model_table.Properties.VariableNames;
    for variable_index = 1:numel(numeric_names)
        variable_name = numeric_names{variable_index};
        values = model_table.(variable_name);
        if iscell(values)
            values = str2double(string(values));
            model_table.(variable_name) = values;
        end
        if isnumeric(values) && ~strcmp(variable_name, x1_name) && ...
                ~strcmp(variable_name, x2_name)
            valid_values = values(~isnan(values));
            unique_values = unique(valid_values);
            if numel(unique_values) > 2 && std(valid_values) > 0
                model_table.(variable_name) = ...
                    (values - mean(valid_values)) ./ std(valid_values);
            end
        end
    end

    predictor_names = {x1_name, x2_name};
    for predictor_index = 1:2
        predictor_name = predictor_names{predictor_index};
        values = model_table.(predictor_name);
        valid_values = values(~isnan(values));
        unique_values = unique(valid_values);
        if numel(unique_values) > 2 && std(valid_values) > 0
            model_table.(predictor_name) = ...
                (values - mean(valid_values)) ./ std(valid_values);
        end
    end

    interaction_name = [x1_name, ':', x2_name];
    model_formula = [y_name, ' ~ ', x1_name, '*', x2_name, ' + ', formula];
    try
        lmm = fitlme(model_table, model_formula);
        coefficient_table = lmm.Coefficients;
        target = coefficient_table(strcmp( ...
            coefficient_table.Name, interaction_name), :);
        if isempty(target)
            reverse_name = [x2_name, ':', x1_name];
            target = coefficient_table(strcmp( ...
                coefficient_table.Name, reverse_name), :);
        end
        if isempty(target)
            error('Interaction coefficient was not found.');
        end
        partial_r2 = target.tStat^2 / (target.tStat^2 + target.DF);
        result_cells{model_index} = {x1_name, x2_name, y_name, ...
            target.tStat, target.pValue, target.DF, target.Estimate, ...
            target.Lower, target.Upper, partial_r2};
    catch model_error
        warning('Interaction model failed for %s: %s', y_name, ...
            model_error.message);
        result_cells{model_index} = {x1_name, x2_name, y_name, ...
            NaN, NaN, NaN, NaN, NaN, NaN, NaN};
    end
end

results = cell2table(vertcat(result_cells{:}), 'VariableNames', ...
    {'x1name', 'x2name', 'yname', 't', 'p', 'DF', 'std_beta', ...
     'ci1', 'ci2', 'Partial_R2'});
end

function variables = extract_formula_variables(formula)
tokens = regexp(formula, '[A-Za-z]\w*', 'match');
keywords = {'fitlme'};
variables = unique(tokens(~ismember(tokens, keywords)), 'stable');
end
