%% config_example.m
% Copy this file to config.m and edit only local paths and run settings.
% Never commit participant-level UK Biobank data or a populated config.m.

PROJECT_ROOT = fileparts(mfilename('fullpath'));
DATA_ROOT = '/path/to/authorized/UKB/data';
OUTPUT_DIR = fullfile(PROJECT_ROOT, 'derived');

TABULAR_DIR = fullfile(DATA_ROOT, 'tabular');
GENETICS_DIR = fullfile(DATA_ROOT, 'genetics');
ATLAS_DIR = fullfile(DATA_ROOT, 'atlas');

FMRI_FILE = fullfile(DATA_ROOT, 'imaging', 'ukb_glasser_ts_2_0.mat');
SC_FILE = fullfile(DATA_ROOT, 'imaging', ...
    'ukb_glasser_s1_tract_2_0_streamline_count.mat');
YEO7_GLASSER_FILE = fullfile(ATLAS_DIR, 'Yeo7MMP.mat');

% Primary parameters reported in the manuscript.
WINDOW_LENGTH = 30;
WINDOW_STEP = 5;
WINDOW_LENGTHS = [30 40 50];
N_WORKERS = 40;
N_CORTICAL_PARCELS = 360;
EXCLUDE_INCOMPLETE_COUPLING = true;

% The current manuscript says that no additional streamline weighting was
% applied. Keep this false to follow the Methods text. The archived analysis
% normalized streamline counts by parcel-pair surface size; set this true only
% after confirming that the archived behavior is the intended paper method.
NORMALIZE_STREAMLINES_BY_PARCEL_SIZE = false;
PARCEL_SIZE_FILE = fullfile(ATLAS_DIR, 'glasser_360_parcel_sizes.mat');

% Analysis-ready participant tables used by scripts 04-08.
BASIC_FILE = fullfile(TABULAR_DIR, 'UKB_Basic.csv');
MRI_COVARIATE_FILE = fullfile(TABULAR_DIR, 'UKB_BrainMRICov.csv');
OUTCOME_FILE = fullfile(TABULAR_DIR, 'UKB_Outcomes_2_0.csv');
LIFESTYLE_FILE = fullfile(TABULAR_DIR, 'UKB_Lifestyle_0_0.csv');
GENETIC_COVARIATE_FILE = fullfile(GENETICS_DIR, 'UKB_GeneCov.csv');
PRS_FILE = fullfile(GENETICS_DIR, 'UKB_PRS.tsv');
APOE_FILE = fullfile(GENETICS_DIR, 'APOE4_genetype.txt');
AD_VARIANT_FILE = fullfile(GENETICS_DIR, 'AD_genetype.txt');

% Optional robustness inputs. Leave empty to skip a comparison.
FA_SFC_FILE = '';
SCHAEFER_SFC_FILE = '';
SCHAEFER_YEO7_FILE = '';

% Optional GWAS settings used by 07_prepare_gwas_inputs.m and 07_run_gwas.sh.
PLINK2_BIN = '/path/to/plink2';
PLINK_BFILE_PATTERN = '/path/to/qc/ukb_imp_chr%d';
GWAS_COVARIATE_FILE = fullfile(GENETICS_DIR, 'UKB_GeneCov_MRI.txt');

if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end
