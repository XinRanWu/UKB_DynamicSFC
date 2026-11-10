%% PRS and APOE associations (Figure 6; Tables S11-S16)
cfg = load_config();
load(cfg.analysis_coupling_table_file,'coupling_360','coupling_7n');
basic = read_analysis_table(cfg.basic_file,'basic/MRI covariates');
gene_cov = read_analysis_table(cfg.genetic_covariate_file,'genetic covariates');
prs = read_analysis_table(cfg.prs_file,'polygenic risk scores');
apoe = read_analysis_table(cfg.apoe_file,'APOE genotype');
ad_variants = read_analysis_table(cfg.ad_variant_file,'candidate AD variants');

pc_names = arrayfun(@(x)sprintf('x22009_0_%d',x),1:10,'UniformOutput',false);
missing_pc = setdiff(pc_names,gene_cov.Properties.VariableNames);
if ~isempty(missing_pc), error('Missing ancestry PCs: %s',strjoin(missing_pc,', ')); end
covariates = [{'AgeAttend_2_0','AgeSq','Sex_0_0','BMI_2_0','HeadMotion_2_0', ...
    'SNR_2_0','Vol_WB_TIV_2_0'},pc_names];

network_data = innerjoin_many(basic,gene_cov,prs,coupling_7n);
network_data(network_data.HeadMotion_2_0>0.20,:) = [];
network_data.AgeSq = zscore(network_data.AgeAttend_2_0).^2;
prs_manifest = readtable(cfg.prs_manifest);
prs_ids = prs_manifest{:,1};
prs_predictors = cell(numel(prs_ids),1);
for i = 1:numel(prs_ids)
    candidates = {sprintf('p%d',prs_ids(i)),sprintf('x%d',prs_ids(i))};
    match = candidates(ismember(candidates,prs.Properties.VariableNames));
    if numel(match) ~= 1
        error('Could not identify one PRS column for UKB field %d.',prs_ids(i));
    end
    prs_predictors{i} = match{1};
end
if numel(prs_predictors) ~= 90
    error('The final PRS manifest must contain 90 scores; found %d.',numel(prs_predictors));
end
network_outcomes = coupling_7n.Properties.VariableNames(2:end);
prs_network = [ ...
    run_lmm_grid(network_data,prs_predictors,network_outcomes(startsWith(network_outcomes,'SFC_')), ...
        covariates,'Centre_0_0',struct('correction','bonferroni','family_size',90*7)); ...
    run_lmm_grid(network_data,prs_predictors,network_outcomes(startsWith(network_outcomes,'DSFC_')), ...
        covariates,'Centre_0_0',struct('correction','bonferroni','family_size',90*7))];
prs_network.correction_method = repmat("Bonferroni within coupling type",height(prs_network),1);
prs_network.correction_family_size = repmat(90*7,height(prs_network),1);

% Regional follow-up of AD-PRS. Configure its exact column name locally if
% the UKB release uses a different label.
ad_candidates = prs_predictors(contains(prs_predictors,'26206','IgnoreCase',true) | ...
    contains(prs_predictors,'alzheimer','IgnoreCase',true));
if numel(ad_candidates) ~= 1
    error('Expected exactly one AD-PRS column; found: %s',strjoin(ad_candidates,', '));
end
parcel_data = innerjoin_many(basic,gene_cov,prs(:,[{'eid'},ad_candidates(:)']),coupling_360);
parcel_data(parcel_data.HeadMotion_2_0>0.20,:) = [];
parcel_data.AgeSq = zscore(parcel_data.AgeAttend_2_0).^2;
parcel_outcomes = coupling_360.Properties.VariableNames(2:end);
ad_prs_parcel = [ ...
    run_lmm_grid(parcel_data,ad_candidates,parcel_outcomes(startsWith(parcel_outcomes,'SFC_')), ...
        covariates,'Centre_0_0',struct('correction','bonferroni','family_size',360)); ...
    run_lmm_grid(parcel_data,ad_candidates,parcel_outcomes(startsWith(parcel_outcomes,'DSFC_')), ...
        covariates,'Centre_0_0',struct('correction','bonferroni','family_size',360))];

variant_predictors = ad_variants.Properties.VariableNames(2:end);
if numel(variant_predictors) ~= 15
    error('The final candidate-variant table must contain 15 SNP dosages; found %d.',numel(variant_predictors));
end
variant_data = innerjoin_many(basic,gene_cov,ad_variants,coupling_7n);
variant_data(variant_data.HeadMotion_2_0>0.20,:) = [];
variant_data.AgeSq = zscore(variant_data.AgeAttend_2_0).^2;
variant_network = [ ...
    run_lmm_grid(variant_data,variant_predictors,network_outcomes(startsWith(network_outcomes,'SFC_')), ...
        covariates,'Centre_0_0',struct('correction','bonferroni','family_size',15*7)); ...
    run_lmm_grid(variant_data,variant_predictors,network_outcomes(startsWith(network_outcomes,'DSFC_')), ...
        covariates,'Centre_0_0',struct('correction','bonferroni','family_size',15*7))];

apoe_network_data = innerjoin_many(basic,gene_cov,apoe,coupling_7n);
apoe_network_data(apoe_network_data.HeadMotion_2_0>0.20,:) = [];
apoe_network_data.AgeSq = zscore(apoe_network_data.AgeAttend_2_0).^2;
apoe_network = [ ...
    run_lmm_grid(apoe_network_data,{'APOE4_dosage'},network_outcomes(startsWith(network_outcomes,'SFC_')), ...
        covariates,'Centre_0_0',struct('correction','bonferroni','family_size',7)); ...
    run_lmm_grid(apoe_network_data,{'APOE4_dosage'},network_outcomes(startsWith(network_outcomes,'DSFC_')), ...
        covariates,'Centre_0_0',struct('correction','bonferroni','family_size',7))];

apoe_data = innerjoin_many(basic,gene_cov,apoe,coupling_360);
apoe_data(apoe_data.HeadMotion_2_0>0.20,:) = [];
apoe_data.AgeSq = zscore(apoe_data.AgeAttend_2_0).^2;
apoe_parcel = [ ...
    run_lmm_grid(apoe_data,{'APOE4_dosage'},parcel_outcomes(startsWith(parcel_outcomes,'SFC_')), ...
        covariates,'Centre_0_0',struct('correction','bonferroni','family_size',360)); ...
    run_lmm_grid(apoe_data,{'APOE4_dosage'},parcel_outcomes(startsWith(parcel_outcomes,'DSFC_')), ...
        covariates,'Centre_0_0',struct('correction','bonferroni','family_size',360))];

output_dir = fullfile(cfg.output_root,'results','06_genetics');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
writetable(prs_network,fullfile(output_dir,'prs_network.csv'));
writetable(ad_prs_parcel,fullfile(output_dir,'ad_prs_parcel.csv'));
writetable(variant_network,fullfile(output_dir,'ad_candidate_variants_network.csv'));
writetable(apoe_network,fullfile(output_dir,'apoe_network.csv'));
writetable(apoe_parcel,fullfile(output_dir,'apoe_parcel.csv'));
