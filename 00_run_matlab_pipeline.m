%% 00_run_matlab_pipeline.m
% Run the primary MATLAB workflow in manuscript order.
% R-based extraction and mediation are run separately from the shell.

if ~exist('config.m', 'file')
    error('Copy config_example.m to config.m and configure local paths.');
end
run('config.m');

fprintf('Stage 2: computing primary %d-TR SFC/DSFC.\n', WINDOW_LENGTH);
run('02_compute_sfc_dsfc_glasser.m');

fprintf('Stage 3: preparing analysis tables.\n');
run('03_prepare_sfc_analysis_tables.m');

fprintf('Stage 4: demographic and age analyses.\n');
run('04_run_demographic_age_associations.m');

fprintf('Stage 5: phenotype and lifestyle analyses.\n');
run('05_run_phenotype_lifestyle_associations.m');

fprintf('Stage 6: PRS, APOE, and AD-variant analyses.\n');
run('06_run_genetic_associations.m');

fprintf('Primary MATLAB workflow completed.\n');
