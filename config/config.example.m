%% Local configuration for the UKB Dynamic SFC analysis
% Copy this file to config/config.m and edit only paths and run settings.
% config/config.m and all derived participant-level data are ignored by Git.

cfg.project_root = fileparts(fileparts(mfilename('fullpath')));
cfg.data_root = '/path/to/authorized/UKB/data';
cfg.output_root = fullfile(cfg.project_root, 'derived');
cfg.analysis_coupling_table_file = fullfile(cfg.output_root, 'tables', 'coupling_tables.mat');

% Imaging inputs. Each MAT file must contain the variables documented in
% docs/DATA_CONTRACT.md. Participant identifiers are matched by eid.
cfg.fmri_file = fullfile(cfg.data_root, 'imaging', 'glasser360_timeseries.mat');
cfg.fmri_id_variable = 'fMRI_Glasser_eID_2_0';
cfg.fmri_data_variable = 'fMRI_Glasser_TS_2_0';
cfg.sc_streamline_file = fullfile(cfg.data_root, 'imaging', 'glasser360_streamline_count.mat');
cfg.sc_fa_file = fullfile(cfg.data_root, 'imaging', 'glasser360_mean_fa.mat');
cfg.sc_id_variable = 'dMRI_Glasser_S1_eID_2_0';
cfg.sc_streamline_variable = 'connectome_streamline_count_10M_2_0';
cfg.sc_fa_variable = 'connectome_mean_fa_2_0';
cfg.gsr_fmri_file = fullfile(cfg.data_root, 'imaging', 'glasser360_timeseries_gsr.mat');

cfg.yeo7_glasser_file = fullfile(cfg.data_root, 'atlas', 'Yeo7MMP.mat');
cfg.yeo7_variable = 'Yeo7MMP';
cfg.spin_indices_file = fullfile(cfg.data_root, 'atlas', 'glasser360_spin_indices_10000.mat');
cfg.spin_indices_variable = 'spin_indices';
cfg.n_parcels = 360;
cfg.network_names = {'VN','SMN','DAN','VAN','LN','FPN','DMN'};

% Primary coupling definition reported in the manuscript.
cfg.window_length = 30;
cfg.window_step = 5;
cfg.edge_rule = 'all'; % all | positive_fc | nonzero_sc | both
cfg.n_workers = 40;
cfg.exclude_incomplete_coupling = true;

% Analysis-ready tables. See docs/DATA_CONTRACT.md.
cfg.basic_file = fullfile(cfg.data_root, 'tabular', 'UKB_BasicMRI.csv');
cfg.outcome_file = fullfile(cfg.data_root, 'tabular', 'UKB_Outcomes_2_0.csv');
cfg.lifestyle_file = fullfile(cfg.data_root, 'tabular', 'UKB_Lifestyle_0_0.csv');
cfg.ses_file = fullfile(cfg.data_root, 'tabular', 'UKB_SES_encoded_0.csv');
cfg.phenotype_manifest = fullfile(cfg.project_root, 'metadata', 'phenotypes_used.csv');
cfg.lifestyle_manifest = fullfile(cfg.project_root, 'metadata', 'lifestyle_used.csv');
cfg.prs_file = fullfile(cfg.data_root, 'genetics', 'UKB_PRS.tsv');
cfg.prs_manifest = fullfile(cfg.project_root, 'UKB_PRS_info.xlsx');
cfg.genetic_covariate_file = fullfile(cfg.data_root, 'genetics', 'UKB_GeneCov.csv');
cfg.apoe_file = fullfile(cfg.data_root, 'genetics', 'APOE4_genotype.csv');
cfg.ad_variant_file = fullfile(cfg.data_root, 'genetics', 'AD_candidate_variants.csv');

% Optional Schaefer-200 sensitivity inputs. Leave empty to skip.
cfg.schaefer_fmri_file = '';
cfg.schaefer_sc_file = '';
cfg.schaefer_yeo7_file = '';
cfg.schaefer_n_parcels = 200;
cfg.schaefer_fmri_id_variable = 'fMRI_Schaefer_eid';
cfg.schaefer_fmri_data_variable = 'fMRI_Schaefer_TS';
cfg.schaefer_sc_id_variable = 'dMRI_Schaefer_eid';
cfg.schaefer_sc_variable = 'Schaefer_streamline_count';
cfg.schaefer_yeo7_variable = 'Yeo7Schaefer';
% Optional parcel-to-fsaverage label vectors for the cross-atlas group map.
cfg.glasser_fsaverage_labels_file = '';
cfg.glasser_fsaverage_labels_variable = 'parcel_labels';
cfg.schaefer_fsaverage_labels_file = '';
cfg.schaefer_fsaverage_labels_variable = 'parcel_labels';
cfg.fsaverage_spin_indices_file = '';
cfg.fsaverage_spin_indices_variable = 'spin_indices';

% GWAS inputs. The source genetic-covariate table must contain UKB Field
% 22021 and Fields 22009.1-22009.20 in addition to eid. The preparation
% script writes the final phenotype, covariate, and --keep files below
% derived/gwas; PLINK applies variant QC within that retained sample.
cfg.plink2_bin = '/path/to/plink2';
cfg.gwas_bfile_pattern = '/path/to/qc/ukb_imp_chr%d';
cfg.gwas_relatedness_variable = 'x22021_0_0';
cfg.gwas_centre_variable = 'Centre_2_0';
cfg.gwas_pc_variables = arrayfun(@(x)sprintf('x22009_0_%d',x),1:20, ...
    'UniformOutput',false);

if ~exist(cfg.output_root, 'dir')
    mkdir(cfg.output_root);
end
