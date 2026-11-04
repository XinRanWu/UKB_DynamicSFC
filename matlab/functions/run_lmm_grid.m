function results = run_lmm_grid(data, predictors, outcomes, covariates, random_intercept, options)
%RUN_LMM_GRID Fit one mixed model for every predictor-outcome pair.
%
% Outcomes are z-scored. Continuous predictors and continuous covariates are
% z-scored within each model's complete-case sample. Binary 0/1 predictors
% retain their original scale. Categorical covariates must be supplied as
% categorical table variables before calling this function.

arguments
    data table
    predictors cell
    outcomes cell
    covariates cell
    random_intercept (1,:) char
    options.correction (1,:) char = 'none'
    options.family_size (1,1) double = NaN
    options.use_parallel (1,1) logical = true
end

required = unique([predictors(:); outcomes(:); covariates(:); {random_intercept}]);
missing = setdiff(required, data.Properties.VariableNames);
if ~isempty(missing)
    error('Missing model variables: %s', strjoin(missing, ', '));
end

n_models = numel(predictors) * numel(outcomes);
rows = cell(n_models,1);
if options.use_parallel
    parfor model_index = 1:n_models
        rows{model_index} = fit_one(model_index, data, predictors, outcomes, ...
            covariates, random_intercept);
    end
else
    for model_index = 1:n_models
        rows{model_index} = fit_one(model_index, data, predictors, outcomes, ...
            covariates, random_intercept);
    end
end
results = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'predictor','outcome','N','beta','SE','t','DF','p','CI_lower','CI_upper','partial_R2'});

switch lower(options.correction)
    case 'bonferroni'
        family_size = options.family_size;
        if ~isfinite(family_size)
            family_size = height(results);
        end
        results.p_corrected = min(results.p * family_size, 1);
    case 'fdr'
        results.p_corrected = mafdr(results.p, 'BHFDR', true);
    case 'none'
        results.p_corrected = results.p;
    otherwise
        error('Unknown correction: %s', options.correction);
end
end

function row = fit_one(model_index, data, predictors, outcomes, covariates, random_intercept)
[outcome_index, predictor_index] = ind2sub([numel(outcomes),numel(predictors)], model_index);
x = predictors{predictor_index};
y = outcomes{outcome_index};
names = unique([{x,y},covariates,{random_intercept}], 'stable');
tb = rmmissing(data(:,names));
try
    tb.(y) = safe_zscore(tb.(y), y);
    x_values = tb.(x);
    is_binary = isnumeric(x_values) && all(ismember(unique(x_values),[0 1]));
    if isnumeric(x_values) && ~is_binary
        tb.(x) = safe_zscore(x_values, x);
    end
    for i = 1:numel(covariates)
        name = covariates{i};
        if isnumeric(tb.(name))
            observed = unique(tb.(name));
            if numel(observed) > 2
                tb.(name) = safe_zscore(tb.(name), name);
            end
        end
    end
    fixed_terms = strjoin([{x},covariates], ' + ');
    formula = sprintf('%s ~ %s + (1|%s)', y, fixed_terms, random_intercept);
    fit = fitlme(tb, formula, 'FitMethod', 'REML');
    coefficient = fit.Coefficients(strcmp(fit.Coefficients.Name,x),:);
    if height(coefficient) ~= 1
        error('Target coefficient %s was not uniquely estimated.', x);
    end
    partial_r2 = coefficient.tStat^2 / (coefficient.tStat^2 + coefficient.DF);
    row = {x,y,height(tb),coefficient.Estimate,coefficient.SE, ...
        coefficient.tStat,coefficient.DF,coefficient.pValue, ...
        coefficient.Lower,coefficient.Upper,partial_r2};
catch model_error
    warning('Model %s ~ %s failed: %s', y, x, model_error.message);
    row = {x,y,height(tb),NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN};
end
end

function z = safe_zscore(x, name)
scale = std(x,0,'omitnan');
if ~isfinite(scale) || scale == 0
    error('%s has zero or invalid variance.', name);
end
z = (x-mean(x,'omitnan'))/scale;
end
