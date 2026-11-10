%% Prepare the final imaging-genetics sample and PLINK 2 input files
% This script implements the sample and covariate definitions reported in
% the revised manuscript. Variant-level MAF, missingness, and HWE filters are
% applied by shell/07_run_gwas.sh after --keep restricts each chromosome to
% this unrelated imaging sample.

cfg = load_config();
load(cfg.analysis_coupling_table_file,'coupling_7n');
basic = read_analysis_table(cfg.basic_file,'basic/MRI covariates');
gene_cov = read_analysis_table(cfg.genetic_covariate_file,'genetic covariates');

phenotypes = [compose('SFC_%s',cfg.network_names), ...
    compose('DSFC_%s',cfg.network_names)];
phenotypes = cellstr(phenotypes(:))';
pc_source = cfg.gwas_pc_variables;
relatedness = cfg.gwas_relatedness_variable;
centre = cfg.gwas_centre_variable;

required_basic = {'eid','Sex_0_0','AgeAttend_2_0','BMI_2_0', ...
    'HeadMotion_2_0','SNR_2_0','Vol_WB_TIV_2_0',centre};
required_gene = [{'eid',relatedness},pc_source];
missing_basic = setdiff(required_basic,basic.Properties.VariableNames);
missing_gene = setdiff(required_gene,gene_cov.Properties.VariableNames);
missing_phenotypes = setdiff(phenotypes,coupling_7n.Properties.VariableNames);
if ~isempty(missing_basic)
    error('Missing GWAS covariates in the basic table: %s',strjoin(missing_basic,', '));
end
if ~isempty(missing_gene)
    error('Missing GWAS genetic variables: %s',strjoin(missing_gene,', '));
end
if ~isempty(missing_phenotypes)
    error('Missing GWAS phenotypes: %s',strjoin(missing_phenotypes,', '));
end

gwas = innerjoin_many(basic(:,required_basic),gene_cov(:,required_gene), ...
    coupling_7n(:,[{'eid'},phenotypes]));

% Field 22021 = 0 indicates that the participant was not identified as
% genetically related to another UK Biobank participant.
gwas = gwas(gwas.(relatedness)==0,:);
gwas = gwas(gwas.HeadMotion_2_0<=0.20,:);
model_variables = [required_basic(2:end),pc_source,phenotypes];
gwas = rmmissing(gwas,'DataVariables',model_variables);
if isempty(gwas)
    error('No complete unrelated imaging participants remained for GWAS.');
end

gwas.FID = gwas.eid;
gwas.IID = gwas.eid;
pheno = gwas(:,[{'FID','IID'},phenotypes]);

age_centered = gwas.AgeAttend_2_0 - mean(gwas.AgeAttend_2_0,'omitnan');
covar = gwas(:,{'FID','IID','Sex_0_0','AgeAttend_2_0','BMI_2_0', ...
    'HeadMotion_2_0','SNR_2_0','Vol_WB_TIV_2_0'});
covar.AgeAttend_2_0_c2 = age_centered.^2;
covar = movevars(covar,'AgeAttend_2_0_c2','After','AgeAttend_2_0');

% PLINK 2 treats nonnumeric strings as categorical covariates and creates
% the required indicator variables internally. Prefixing the centre code
% prevents it from being interpreted as a continuous numeric covariate.
covar.ImagingCentre = categorical("site_"+string(gwas.(centre)));
for i = 1:numel(pc_source)
    covar.(sprintf('PC%d',i)) = gwas.(pc_source{i});
end

keep = gwas(:,{'FID','IID'});
output_dir = fullfile(cfg.output_root,'gwas');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
writetable(pheno,fullfile(output_dir,'coupling_phenotypes.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(covar,fullfile(output_dir,'gwas_covariates.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(keep,fullfile(output_dir,'unrelated_imaging_keep.tsv'), ...
    'FileType','text','Delimiter','\t','WriteVariableNames',false);

sample_summary = table(height(gwas),mean(gwas.AgeAttend_2_0), ...
    std(gwas.AgeAttend_2_0),sum(gwas.Sex_0_0==0),sum(gwas.Sex_0_0==1), ...
    'VariableNames',{'N','Age_mean','Age_SD','Sex_0_N','Sex_1_N'});
writetable(sample_summary,fullfile(output_dir,'gwas_sample_summary.csv'));
