%% 07_prepare_gwas_inputs.m
% Create the PLINK phenotype table for network-level GWAS.

if ~exist('config.m', 'file')
    error('Copy config_example.m to config.m and configure local paths.');
end
run('config.m');
window_length = WINDOW_LENGTH;

brain_file = fullfile(OUTPUT_DIR, 'analysis_tables', sprintf( ...
    'brain_data_%dTR.mat', window_length));
if ~exist(brain_file, 'file')
    error('Run 03_prepare_sfc_analysis_tables.m first.');
end
load(brain_file, 'brain_data_7n', 'brain_data_360p');

phenotype_table = brain_data_7n;
phenotype_table.SFC_Global = mean( ...
    brain_data_360p{:, startsWith(brain_data_360p.Properties.VariableNames, ...
    'SFC_P')}, 2, 'omitnan');
phenotype_table.DSFC_Global = mean( ...
    brain_data_360p{:, startsWith(brain_data_360p.Properties.VariableNames, ...
    'DSFC_P')}, 2, 'omitnan');

phenotype_table = movevars(phenotype_table, ...
    'SFC_Global', 'After', 'SFC_Default');
phenotype_table = movevars(phenotype_table, ...
    'DSFC_Global', 'After', 'DSFC_Default');
phenotype_table.FID = phenotype_table.eid;
phenotype_table.IID = phenotype_table.eid;
phenotype_table = movevars(phenotype_table, {'FID', 'IID'}, 'Before', 'eid');
phenotype_table = removevars(phenotype_table, 'eid');

gwas_dir = fullfile(OUTPUT_DIR, 'gwas');
if ~exist(gwas_dir, 'dir')
    mkdir(gwas_dir);
end
phenotype_file = fullfile(gwas_dir, sprintf( ...
    'UKB_sfc_dsfc_phenotypes_%dTR.tsv', window_length));
writetable(phenotype_table, phenotype_file, ...
    'FileType', 'text', 'Delimiter', '\t');
fprintf('GWAS phenotype table saved to %s\n', phenotype_file);
