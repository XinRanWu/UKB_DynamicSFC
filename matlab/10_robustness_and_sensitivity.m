%% Methodological sensitivity analyses (Supplementary Figures S1-S8)
cfg = load_config();
primary_file = fullfile(cfg.output_root,'coupling','primary_glasser360_30TR.mat');
load(primary_file,'primary');
fmri = load_mat_variable(cfg.fmri_file,cfg.fmri_data_variable);
fmri_eid = load_mat_variable(cfg.fmri_file,cfg.fmri_id_variable);
sc = load_mat_variable(cfg.sc_streamline_file,cfg.sc_streamline_variable);
sc_eid = load_mat_variable(cfg.sc_streamline_file,cfg.sc_id_variable);
yeo = load_mat_variable(cfg.yeo7_glasser_file,cfg.yeo7_variable);
spin_indices = [];
if exist(cfg.spin_indices_file,'file')
    spin_indices = load_mat_variable(cfg.spin_indices_file,cfg.spin_indices_variable);
else
    warning('Spin indices not found; group-map correlations will be exported without spin-test P values.');
end
base_options = struct('n_parcels',cfg.n_parcels,'window_length',cfg.window_length, ...
    'window_step',cfg.window_step,'edge_rule','all','n_workers',cfg.n_workers, ...
    'exclude_incomplete',cfg.exclude_incomplete_coupling);

output_dir = fullfile(cfg.output_root,'sensitivity');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
comparisons = {}; network_comparisons = {}; alternatives = struct();

% Alternative sliding-window lengths.
for window = [40 50]
    options = base_options; options.window_length = window;
    name = sprintf('window_%dTR',window);
    alternatives.(name) = compute_sfc_dsfc(fmri,fmri_eid,sc,sc_eid,options);
    [comparisons{end+1},network_comparisons{end+1}] = ...
        compare_coupling_estimates(primary,alternatives.(name),yeo,name,spin_indices); %#ok<SAGROW>
end

% Edge-inclusion definitions.
for rule = {'positive_fc','nonzero_sc','both'}
    options = base_options; options.edge_rule = rule{1};
    name = ['edge_' rule{1}];
    alternatives.(name) = compute_sfc_dsfc(fmri,fmri_eid,sc,sc_eid,options);
    [comparisons{end+1},network_comparisons{end+1}] = ...
        compare_coupling_estimates(primary,alternatives.(name),yeo,name,spin_indices); %#ok<SAGROW>
end

% Fractional-anisotropy structural weights.
if exist(cfg.sc_fa_file,'file')
    fa = load_mat_variable(cfg.sc_fa_file,cfg.sc_fa_variable);
    fa_eid = load_mat_variable(cfg.sc_fa_file,cfg.sc_id_variable);
    alternatives.fa = compute_sfc_dsfc(fmri,fmri_eid,fa,fa_eid,base_options);
    [comparisons{end+1},network_comparisons{end+1}] = ...
        compare_coupling_estimates(primary,alternatives.fa,yeo,'streamline_vs_FA',spin_indices);
end

% Global-signal-regressed fMRI input.
if exist(cfg.gsr_fmri_file,'file')
    gsr = load_mat_variable(cfg.gsr_fmri_file,cfg.fmri_data_variable);
    gsr_eid = load_mat_variable(cfg.gsr_fmri_file,cfg.fmri_id_variable);
    alternatives.gsr = compute_sfc_dsfc(gsr,gsr_eid,sc,sc_eid,base_options);
    [comparisons{end+1},network_comparisons{end+1}] = ...
        compare_coupling_estimates(primary,alternatives.gsr,yeo,'without_vs_with_GSR',spin_indices);
end

% Schaefer-200 parcellation. Parcelwise correspondence is not computed because
% the atlases have no one-to-one parcels. We export matched Yeo-7 correlations
% and, when projection labels are configured, fsaverage group-map spin tests.
if ~isempty(cfg.schaefer_fmri_file) && exist(cfg.schaefer_fmri_file,'file') && ...
        exist(cfg.schaefer_sc_file,'file') && exist(cfg.schaefer_yeo7_file,'file')
    sfmri = load_mat_variable(cfg.schaefer_fmri_file,cfg.schaefer_fmri_data_variable);
    sfmri_eid = load_mat_variable(cfg.schaefer_fmri_file,cfg.schaefer_fmri_id_variable);
    ssc = load_mat_variable(cfg.schaefer_sc_file,cfg.schaefer_sc_variable);
    ssc_eid = load_mat_variable(cfg.schaefer_sc_file,cfg.schaefer_sc_id_variable);
    soptions = base_options; soptions.n_parcels = cfg.schaefer_n_parcels;
    alternatives.schaefer200 = compute_sfc_dsfc(sfmri,sfmri_eid,ssc,ssc_eid,soptions);
    syeo = load_mat_variable(cfg.schaefer_yeo7_file,cfg.schaefer_yeo7_variable);
    [eid,ip,is] = intersect(primary.eid,alternatives.schaefer200.eid,'stable');
    measures = {'SFC','DSFC'};
    for m = 1:2
        p7 = aggregate_yeo7(primary.(measures{m})(ip,:),yeo);
        s7 = aggregate_yeo7(alternatives.schaefer200.(measures{m})(is,:),syeo);
        for network = 1:7
            [r,p] = corr(p7(:,network),s7(:,network),'Rows','pairwise');
            network_comparisons{end+1} = cell2table( ...
                {'glasser360_vs_schaefer200',measures{m},network,r,p}, ...
                'VariableNames',{'comparison','measure','network','pearson_r','p'}); %#ok<SAGROW>
        end
    end
    if exist(cfg.glasser_fsaverage_labels_file,'file') && ...
            exist(cfg.schaefer_fsaverage_labels_file,'file') && ...
            exist(cfg.fsaverage_spin_indices_file,'file')
        glabels = load_mat_variable(cfg.glasser_fsaverage_labels_file,cfg.glasser_fsaverage_labels_variable);
        slabels = load_mat_variable(cfg.schaefer_fsaverage_labels_file,cfg.schaefer_fsaverage_labels_variable);
        fspin = load_mat_variable(cfg.fsaverage_spin_indices_file,cfg.fsaverage_spin_indices_variable);
        for m = 1:2
            gmap = project_parcel_map(mean(primary.(measures{m}),1,'omitnan'),glabels);
            smap = project_parcel_map(mean(alternatives.schaefer200.(measures{m}),1,'omitnan'),slabels);
            [rho,p_spin] = spin_test_correlation(gmap,smap,fspin);
            comparisons{end+1} = cell2table( ...
                {'glasser360_vs_schaefer200',measures{m},numel(eid),NaN,NaN,rho,p_spin,'all'}, ...
                'VariableNames',{'comparison','measure','N','mean_within_spearman', ...
                'median_within_spearman','group_spearman','group_spin_p','edge_rule'}); %#ok<SAGROW>
        end
    end
end

% Intrinsic BOLD properties used as additional covariates downstream.
bold_properties = compute_bold_properties(fmri,fmri_eid,cfg.n_parcels);
save(fullfile(output_dir,'alternative_coupling_and_bold_properties.mat'), ...
    'alternatives','bold_properties','-v7.3');
writetable(vertcat(comparisons{:}),fullfile(output_dir,'individual_group_correspondence.csv'));
writetable(vertcat(network_comparisons{:}),fullfile(output_dir,'network_correspondence.csv'));

% Export every alternative using the same table contract as the primary
% analysis. Set cfg.analysis_coupling_table_file to one of these MAT files
% and rerun scripts 03-06 to reproduce association-level sensitivity tests.
alternative_names = fieldnames(alternatives);
for i = 1:numel(alternative_names)
    name = alternative_names{i};
    labels = yeo;
    if strcmp(name,'schaefer200'), labels = syeo; end
    [coupling_360,coupling_7n] = coupling_to_tables( ...
        alternatives.(name),labels,cfg.network_names);
    variant_dir = fullfile(output_dir,'tables',name);
    if ~exist(variant_dir,'dir'), mkdir(variant_dir); end
    save(fullfile(variant_dir,'coupling_tables.mat'),'coupling_360','coupling_7n','-v7.3');
    writetable(coupling_360,fullfile(variant_dir,'coupling_360.csv'));
    writetable(coupling_7n,fullfile(variant_dir,'coupling_7n.csv'));
end

fprintf(['Alternative coupling estimates are saved for rerunning scripts 02-06. ', ...
    'Use the same scripts after replacing the primary coupling input, so model ', ...
    'specifications and multiple-comparison families remain unchanged.\n']);
