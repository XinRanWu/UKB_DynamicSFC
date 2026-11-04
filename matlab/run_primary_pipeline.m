%% Run primary analyses in the same order as the manuscript Results
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo_root,'matlab','functions'));

run(fullfile(repo_root,'matlab','01_compute_primary_coupling.m'));
run(fullfile(repo_root,'matlab','02_prepare_analysis_tables.m'));
run(fullfile(repo_root,'matlab','03_descriptive_and_age_associations.m'));
run(fullfile(repo_root,'matlab','04_phenotype_associations.m'));
run(fullfile(repo_root,'matlab','05_lifestyle_associations.m'));
run(fullfile(repo_root,'matlab','06_genetic_associations.m'));
run(fullfile(repo_root,'matlab','07_prepare_gwas_inputs.m'));

fprintf('\nPrimary MATLAB pipeline completed. Run scripts 08-11 as documented in README.\n');
