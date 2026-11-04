function results = run_lmm_matched_adjustment(data,predictors,outcomes,base_covariates, ...
    random_intercept,adjustment_prefix)
%RUN_LMM_MATCHED_ADJUSTMENT Add the regional/network-matched BOLD covariate.
results = table();
for i = 1:numel(outcomes)
    outcome = outcomes{i};
    suffix = regexprep(outcome,'^(SFC|DSFC)_','');
    matched_covariate = [adjustment_prefix '_' suffix];
    if ~ismember(matched_covariate,data.Properties.VariableNames)
        error('Missing matched covariate %s for %s.',matched_covariate,outcome);
    end
    piece = run_lmm_grid(data,predictors,{outcome}, ...
        [base_covariates,{matched_covariate}],random_intercept, ...
        struct('correction','none'));
    piece.adjustment(:) = {adjustment_prefix};
    results = [results;piece]; %#ok<AGROW>
end
end
