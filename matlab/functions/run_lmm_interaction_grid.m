function results = run_lmm_interaction_grid(data, predictors, moderator, outcomes, covariates, random_intercept)
%RUN_LMM_INTERACTION_GRID Test predictor by moderator interactions.
n_models = numel(predictors)*numel(outcomes); rows = cell(n_models,1);
parfor model_index = 1:n_models
    [outcome_index,predictor_index] = ind2sub([numel(outcomes),numel(predictors)],model_index);
    x = predictors{predictor_index}; y = outcomes{outcome_index};
    names = unique([{x,moderator,y},covariates,{random_intercept}],'stable');
    tb = rmmissing(data(:,names));
    try
        tb.(y) = z(tb.(y)); tb.(x) = z_unless_binary(tb.(x));
        tb.(moderator) = z_unless_binary(tb.(moderator));
        for i = 1:numel(covariates)
            name = covariates{i}; values = tb.(name);
            if isnumeric(values) && numel(unique(values))>2, tb.(name)=z(values); end
        end
        formula = sprintf('%s ~ %s*%s + %s + (1|%s)',y,x,moderator, ...
            strjoin(covariates,' + '),random_intercept);
        fit = fitlme(tb,formula,'FitMethod','REML');
        names_expected = {[x ':' moderator],[moderator ':' x]};
        coefficient = fit.Coefficients(ismember(fit.Coefficients.Name,names_expected),:);
        if height(coefficient)~=1, error('Interaction coefficient not uniquely estimated.'); end
        rows{model_index} = {x,moderator,y,height(tb),coefficient.Estimate,coefficient.SE, ...
            coefficient.tStat,coefficient.DF,coefficient.pValue,coefficient.Lower,coefficient.Upper};
    catch model_error
        warning('Interaction model %s failed: %s',y,model_error.message);
        rows{model_index} = {x,moderator,y,height(tb),NaN,NaN,NaN,NaN,NaN,NaN,NaN};
    end
end
results = cell2table(vertcat(rows{:}),'VariableNames', ...
    {'predictor','moderator','outcome','N','beta','SE','t','DF','p','CI_lower','CI_upper'});
results.p_corrected = min(results.p * height(results),1);
end

function value = z_unless_binary(value)
if isnumeric(value) && all(ismember(unique(value),[0 1])), return; end
value = z(value);
end

function value = z(value)
s = std(value,0,'omitnan');
if ~isfinite(s)||s==0, error('Variable has zero or invalid variance.'); end
value = (value-mean(value,'omitnan'))/s;
end
