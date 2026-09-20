function gwas = prepare_gwas_phenotypes(window_length)
% PREPARE_GWAS_PHENOTYPES
%
% Prepare participant-level SFC/DSFC phenotypes and covariates for PLINK2 GWAS.
%
% This script reproduces the phenotype organization used in the original
% master script while separating phenotype preparation from chromosome-level
% PLINK execution.
%
% Required inputs:
%   brain_data_<window>TR.mat
%   UKB_GWASCov.csv
%
% UKB_GWASCov.csv should contain:
%   eid
%   Sex_0_0
%   AgeAttend_2_0
%   BMI_2_0
%   HeadMotion_2_0
%   Vol_WB_TIV_2_0
%   Centre_2_0  (or Centre_0_0 if that is the verified study variable)
%   PC1 ... PC20
%
% The original code used network-level and global SFC/DSFC phenotypes.
%
% -------------------------------------------------------------------------

if nargin < 1 || isempty(window_length)
    window_length = 30;
end

run('00_config.m');

table_dir = fullfile(OUTPUT_DIR, 'analysis_tables');

brain_file = fullfile( ...
    table_dir, ...
    sprintf('brain_data_%dTR.mat', window_length));

if ~exist(brain_file, 'file')
    error('Run prepare_sfc_analysis_tables.m first.');
end

load(brain_file, ...
    'brain_data_7n', ...
    'brain_data_360p');

% -------------------------------------------------------------------------
% Reconstruct global SFC/DSFC from parcel-level values
% -------------------------------------------------------------------------

sfc_vars = brain_data_360p.Properties.VariableNames( ...
    startsWith(brain_data_360p.Properties.VariableNames, 'SFC_P'));

dsfc_vars = brain_data_360p.Properties.VariableNames( ...
    startsWith(brain_data_360p.Properties.VariableNames, 'DSFC_P'));

SFC_Global = mean( ...
    brain_data_360p{:, sfc_vars}, ...
    2, ...
    'omitnan');

DSFC_Global = mean( ...
    brain_data_360p{:, dsfc_vars}, ...
    2, ...
    'omitnan');

% -------------------------------------------------------------------------
% Build GWAS phenotype table
% -------------------------------------------------------------------------

network_names = { ...
    'Visual', ...
    'Somatomotor', ...
    'DorsalAttention', ...
    'VentralAttention', ...
    'Limbic', ...
    'Frontoparietal', ...
    'Default' ...
};

phenotype = table();

phenotype.FID = brain_data_7n.eid;
phenotype.IID = brain_data_7n.eid;

for n = 1:7
    src = sprintf('SFC_%s', network_names{n});
    dst = sprintf('SFC_%s', gwas_abbreviation(network_names{n}));
    phenotype.(dst) = brain_data_7n.(src);
end

phenotype.SFC_Global = SFC_Global;

for n = 1:7
    src = sprintf('DSFC_%s', network_names{n});
    dst = sprintf('DSFC_%s', gwas_abbreviation(network_names{n}));
    phenotype.(dst) = brain_data_7n.(src);
end

phenotype.DSFC_Global = DSFC_Global;

% -------------------------------------------------------------------------
% Apply head-motion threshold using verified GWAS covariate table
% -------------------------------------------------------------------------

cov_file = fullfile( ...
    DATA_ROOT, ...
    'genetics', ...
    'UKB_GWASCov.csv');

if ~exist(cov_file, 'file')
    error(['GWAS covariate file not found: ', cov_file, newline, ...
           'Create this table from the verified study covariates before ', ...
           'running the GWAS preparation script.']);
end

cov = readtable(cov_file);

cov = harmonize_eid(cov);

if ~ismember('HeadMotion_2_0', cov.Properties.VariableNames)
    error('UKB_GWASCov.csv must contain HeadMotion_2_0.');
end

keep_ids = cov.eid( ...
    isnan(cov.HeadMotion_2_0) | cov.HeadMotion_2_0 <= 0.2);

phenotype = phenotype(ismember(phenotype.IID, keep_ids), :);

% -------------------------------------------------------------------------
% Export PLINK-compatible phenotype file
% -------------------------------------------------------------------------

gwas_dir = fullfile(OUTPUT_DIR, 'gwas');

if ~exist(gwas_dir, 'dir')
    mkdir(gwas_dir);
end

pheno_file = fullfile( ...
    gwas_dir, ...
    sprintf('UKB_sfc_dsfc_pheno_%dTR.txt', window_length));

writetable( ...
    phenotype, ...
    pheno_file, ...
    'Delimiter', ...
    ' ', ...
    'FileType', ...
    'text');

% -------------------------------------------------------------------------
% Export matching covariate table
% -------------------------------------------------------------------------

cov = cov(ismember(cov.eid, phenotype.IID), :);

% Match phenotype ordering.
[tf, order] = ismember(phenotype.IID, cov.eid);

if ~all(tf)
    error('Some phenotype participants are missing from the GWAS covariates.');
end

cov = cov(order, :);

if ismember('FID', cov.Properties.VariableNames)
    cov.FID = phenotype.FID;
else
    cov.FID = phenotype.FID;
end

if ismember('IID', cov.Properties.VariableNames)
    cov.IID = phenotype.IID;
else
    cov.IID = phenotype.IID;
end

% Put FID/IID first.
remaining = setdiff( ...
    cov.Properties.VariableNames, ...
    {'FID', 'IID', 'eid'}, ...
    'stable');

cov_out = cov(:, [{'FID', 'IID'}, remaining]);

covar_file = fullfile( ...
    gwas_dir, ...
    sprintf('UKB_sfc_dsfc_covar_%dTR.txt', window_length));

writetable( ...
    cov_out, ...
    covar_file, ...
    'Delimiter', ...
    ' ', ...
    'FileType', ...
    'text');

% -------------------------------------------------------------------------
% Save metadata needed by the PLINK execution script
% -------------------------------------------------------------------------

gwas = struct();
gwas.phenotype_file = pheno_file;
gwas.covariate_file = covar_file;
gwas.phenotype_names = phenotype.Properties.VariableNames(3:end);
gwas.n_participants = height(phenotype);

save( ...
    fullfile(gwas_dir, ...
    sprintf('gwas_preparation_%dTR.mat', window_length)), ...
    'gwas');

fprintf('GWAS phenotype file: %s\n', pheno_file);
fprintf('GWAS covariate file: %s\n', covar_file);
fprintf('Participants: %d\n', gwas.n_participants);

end


function abb = gwas_abbreviation(network_name)

switch network_name
    case 'Visual'
        abb = 'VN';
    case 'Somatomotor'
        abb = 'SMN';
    case 'DorsalAttention'
        abb = 'DAN';
    case 'VentralAttention'
        abb = 'VAN';
    case 'Limbic'
        abb = 'LN';
    case 'Frontoparietal'
        abb = 'FPN';
    case 'Default'
        abb = 'DMN';
    otherwise
        error('Unknown network name: %s', network_name);
end

end


function T = harmonize_eid(T)

vars = T.Properties.VariableNames;
candidates = {'eid', 'IID', 'iid', 'ID', 'participant_id'};

for i = 1:numel(candidates)
    if ismember(candidates{i}, vars)
        if ~strcmp(candidates{i}, 'eid')
            T.Properties.VariableNames{strcmp(vars, candidates{i})} = 'eid';
        end
        return;
    end
end

error('No participant-ID column could be identified.');
end
