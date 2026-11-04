%% Repeat principal associations after AC1 or BOLD-variability adjustment
cfg = load_config();
load(cfg.analysis_coupling_table_file,'coupling_360','coupling_7n');
load(fullfile(cfg.output_root,'sensitivity','alternative_coupling_and_bold_properties.mat'), ...
    'bold_properties');
yeo = load_mat_variable(cfg.yeo7_glasser_file,cfg.yeo7_variable);
[bold360,bold7] = bold_properties_to_tables(bold_properties,yeo,cfg.network_names);
basic = read_analysis_table(cfg.basic_file,'basic/MRI covariates');
phenotypes = read_analysis_table(cfg.outcome_file,'phenotypes');
prs = read_analysis_table(cfg.prs_file,'PRS');
gene = read_analysis_table(cfg.genetic_covariate_file,'genetic covariates');

base = {'AgeSq','Sex_0_0','Race_1','Race_2','Race_3','Race_4','BMI_2_0', ...
    'HeadMotion_2_0','Vol_WB_TIV_2_0'};
age_data = innerjoin_many(basic,coupling_360,bold360);
age_data(age_data.HeadMotion_2_0>0.20,:) = [];
age_data.AgeSq = zscore(age_data.AgeAttend_2_0).^2;
parcel_outcomes = coupling_360.Properties.VariableNames(2:end);
age_ac1 = run_lmm_matched_adjustment(age_data,{'AgeAttend_2_0'},parcel_outcomes, ...
    base,'Centre_0_0','AC1');
age_bv = run_lmm_matched_adjustment(age_data,{'AgeAttend_2_0'},parcel_outcomes, ...
    base,'Centre_0_0','BV');

manifest = readtable(cfg.phenotype_manifest,'TextType','char');
phenotype_predictors = manifest.variable(:)';
phenotype_data = innerjoin_many(basic,phenotypes,coupling_7n,bold7);
phenotype_data(phenotype_data.HeadMotion_2_0>0.20,:) = [];
phenotype_data.AgeSq = zscore(phenotype_data.AgeAttend_2_0).^2;
phenotype_base = {'AgeAttend_2_0','AgeSq','Sex_0_0','Race_1','Race_2','Race_3','Race_4', ...
    'BMI_2_0','HeadMotion_2_0','SNR_2_0','Vol_WB_TIV_2_0'};
network_outcomes = coupling_7n.Properties.VariableNames(2:end);
phenotype_ac1 = run_lmm_matched_adjustment(phenotype_data,phenotype_predictors, ...
    network_outcomes,phenotype_base,'Centre_0_0','AC1');
phenotype_bv = run_lmm_matched_adjustment(phenotype_data,phenotype_predictors, ...
    network_outcomes,phenotype_base,'Centre_0_0','BV');

ad_candidates = prs.Properties.VariableNames(contains(prs.Properties.VariableNames,'26206') | ...
    contains(prs.Properties.VariableNames,'alzheimer','IgnoreCase',true));
if numel(ad_candidates)~=1, error('Expected exactly one AD-PRS variable.'); end
pc_names = arrayfun(@(x)sprintf('x22009_0_%d',x),1:10,'UniformOutput',false);
prs_data = innerjoin_many(basic,gene,prs(:,[{'eid'},ad_candidates(:)']),coupling_360,bold360);
prs_data(prs_data.HeadMotion_2_0>0.20,:) = [];
prs_data.AgeSq = zscore(prs_data.AgeAttend_2_0).^2;
prs_base = [{'AgeAttend_2_0','AgeSq','Sex_0_0','BMI_2_0','HeadMotion_2_0', ...
    'SNR_2_0','Vol_WB_TIV_2_0'},pc_names];
prs_ac1 = run_lmm_matched_adjustment(prs_data,ad_candidates,parcel_outcomes, ...
    prs_base,'Centre_0_0','AC1');
prs_bv = run_lmm_matched_adjustment(prs_data,ad_candidates,parcel_outcomes, ...
    prs_base,'Centre_0_0','BV');

output = fullfile(cfg.output_root,'sensitivity','intrinsic_bold_adjustment');
if ~exist(output,'dir'), mkdir(output); end
writetable([age_ac1;age_bv],fullfile(output,'age_associations.csv'));
writetable([phenotype_ac1;phenotype_bv],fullfile(output,'phenotype_associations.csv'));
writetable([prs_ac1;prs_bv],fullfile(output,'ad_prs_associations.csv'));
