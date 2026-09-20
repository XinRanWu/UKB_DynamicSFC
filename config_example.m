%% config_example.m
%
% Example configuration for the Dynamic SFC analysis.
%
% Copy this file to:
%
%     config.m
%
% and modify the paths for your own computing environment.
%
% Do NOT commit participant-level UK Biobank data to the public repository.

%% ------------------------------------------------------------------------
% Data directories
% -------------------------------------------------------------------------

DATA_ROOT = '/path/to/authorized/UKB/data';

FMRI_FILE = fullfile( ...
    DATA_ROOT, ...
    'fMRI', ...
    'ukb_glasser_ts_2_0.mat' ...
);

SC_FILE = fullfile( ...
    DATA_ROOT, ...
    'dMRI', ...
    'ukb_glasser_s1_tract_2_0_streamline_count.mat' ...
);

OUTPUT_DIR = fullfile( ...
    DATA_ROOT, ...
    'derived', ...
    'dynamic_sfc' ...
);

%% ------------------------------------------------------------------------
% Analysis parameters
% -------------------------------------------------------------------------

WINDOW_LENGTH = 30;

WINDOW_STEP = 5;

N_WORKERS = 40;

%% ------------------------------------------------------------------------
% Atlas
% -------------------------------------------------------------------------

ATLAS_NAME = 'glasser_360_fsa5';

N_CORTICAL_PARCELS = 360;

%% ------------------------------------------------------------------------
% Create output directory
% -------------------------------------------------------------------------

if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end
