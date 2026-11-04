%% Environmental and lifestyle associations (Figure 4)
cfg = load_config();
load(cfg.analysis_coupling_table_file,'coupling_360','coupling_7n');
basic = read_analysis_table(cfg.basic_file,'basic/MRI covariates');
lifestyle = read_analysis_table(cfg.lifestyle_file,'lifestyle/environment');

network_data = innerjoin_many(basic,lifestyle,coupling_7n);
network_data(network_data.HeadMotion_2_0>0.20,:) = [];
network_data.AgeSq = zscore(network_data.AgeAttend_2_0).^2;
manifest = readtable(cfg.lifestyle_manifest,'TextType','char');
predictors = manifest.variable(:)';
missing_predictors = setdiff(predictors,lifestyle.Properties.VariableNames);
if ~isempty(missing_predictors)
    error('Lifestyle table is missing final analysis variables: %s',strjoin(missing_predictors,', '));
end
outcomes7 = coupling_7n.Properties.VariableNames(2:end);
covariates = {'AgeAttend_2_0','AgeSq','Sex_0_0','Race_1','Race_2','Race_3','Race_4', ...
    'BMI_2_0','HeadMotion_2_0','SNR_2_0','Vol_WB_TIV_2_0'};

network_results = [ ...
    run_lmm_grid(network_data,predictors,outcomes7(startsWith(outcomes7,'SFC_')),covariates, ...
        'Centre_0_0',struct('correction','bonferroni','family_size',numel(predictors)*7)); ...
    run_lmm_grid(network_data,predictors,outcomes7(startsWith(outcomes7,'DSFC_')),covariates, ...
        'Centre_0_0',struct('correction','bonferroni','family_size',numel(predictors)*7))];
follow_up = unique(network_results.predictor(network_results.p_corrected<0.05));
parcel_results = table();
if ~isempty(follow_up)
    parcel_data = innerjoin_many(basic,lifestyle(:,[{'eid'},follow_up(:)']),coupling_360);
    parcel_data(parcel_data.HeadMotion_2_0>0.20,:) = [];
    parcel_data.AgeSq = zscore(parcel_data.AgeAttend_2_0).^2;
    outcomes360 = coupling_360.Properties.VariableNames(2:end);
    parcel_results = [ ...
        run_lmm_grid(parcel_data,follow_up,outcomes360(startsWith(outcomes360,'SFC_')),covariates, ...
            'Centre_0_0',struct('correction','bonferroni','family_size',360)); ...
        run_lmm_grid(parcel_data,follow_up,outcomes360(startsWith(outcomes360,'DSFC_')),covariates, ...
            'Centre_0_0',struct('correction','bonferroni','family_size',360))];
end

output_dir = fullfile(cfg.output_root,'results','05_lifestyle');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
writetable(network_results,fullfile(output_dir,'lifestyle_network.csv'));
writetable(parcel_results,fullfile(output_dir,'lifestyle_parcel_followup.csv'));
