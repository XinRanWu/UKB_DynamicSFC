function outputs = prepare_sfc_analysis_tables(window_lengths)
% PREPARE_SFC_ANALYSIS_TABLES
%
% Convert participant-level Glasser-360 SFC/DSFC outputs into analysis tables
% used by downstream association analyses.
%
% This script reproduces the organization used in the original master script:
%   - parcel-level SFC/DSFC (360 parcels)
%   - Yeo-7 network-averaged SFC/DSFC
%   - participant ID stored as variable 'eid'
%
% INPUT
% -----
% window_lengths : vector of window lengths in TR.
%                  Default: [30 40 50]
%
% REQUIRED CONFIG
% ---------------
% 00_config.m must define:
%   OUTPUT_DIR
%
% REQUIRED FILE
% -------------
% Yeo7MMP must be available in:
%   atlas/Yeo7MMP.mat
%
% The file should contain a 360 x 1 network-label vector named Yeo7MMP.
%
% OUTPUT
% ------
% A struct with one entry per window length:
%   outputs.TR30.brain_data_7n
%   outputs.TR30.brain_data_360p
%   ...
%
% In addition, analysis tables are saved to derived/analysis_tables.
%
% -------------------------------------------------------------------------

if nargin < 1 || isempty(window_lengths)
    window_lengths = [30 40 50];
end

if ~exist('00_config.m', 'file')
    error(['00_config.m was not found. Copy 00_config_example.m to ', ...
           '00_config.m and configure the paths.']);
end

run('00_config.m');

analysis_table_dir = fullfile(OUTPUT_DIR, 'analysis_tables');

if ~exist(analysis_table_dir, 'dir')
    mkdir(analysis_table_dir);
end

atlas_file = fullfile('atlas', 'Yeo7MMP.mat');

if ~exist(atlas_file, 'file')
    error(['Yeo7MMP atlas mapping not found: ', atlas_file]);
end

load(atlas_file, 'Yeo7MMP');

if numel(Yeo7MMP) ~= 360
    error('Yeo7MMP must contain one Yeo-7 network label for each Glasser parcel.');
end

Yeo7MMP = Yeo7MMP(:);

Yeo7Names = { ...
    'Visual'; ...
    'Somatomotor'; ...
    'DorsalAttention'; ...
    'VentralAttention'; ...
    'Limbic'; ...
    'Frontoparietal'; ...
    'Default' ...
};

outputs = struct();

for w = window_lengths

    input_file = fullfile( ...
        OUTPUT_DIR, ...
        sprintf('ukb_glasser_scfcc_2_0_streamline_count_%dTR.mat', w));

    if ~exist(input_file, 'file')
        warning('Input not found and will be skipped: %s', input_file);
        continue;
    end

    load(input_file, ...
        'SCDFCC_Mean', ...
        'SCDFCC_SD', ...
        'common_eid');

    if size(SCDFCC_Mean, 2) ~= 360 || size(SCDFCC_SD, 2) ~= 360
        error('Expected 360 Glasser parcels in %s.', input_file);
    end

    if size(SCDFCC_Mean, 1) ~= numel(common_eid) || ...
       size(SCDFCC_SD, 1) ~= numel(common_eid)
        error('Participant count does not match common_eid in %s.', input_file);
    end

    % ---------------------------------------------------------------------
    % Yeo-7 network means
    % ---------------------------------------------------------------------

    SFC_7n = network_mean(SCDFCC_Mean, Yeo7MMP, 7);
    DSFC_7n = network_mean(SCDFCC_SD, Yeo7MMP, 7);

    brain_data_7n = table(common_eid(:), 'VariableNames', {'eid'});

    for n = 1:7
        brain_data_7n.(sprintf('SFC_%s', Yeo7Names{n})) = SFC_7n(:, n);
    end

    for n = 1:7
        brain_data_7n.(sprintf('DSFC_%s', Yeo7Names{n})) = DSFC_7n(:, n);
    end

    % ---------------------------------------------------------------------
    % Parcel-level table
    % ---------------------------------------------------------------------

    brain_data_360p = table(common_eid(:), 'VariableNames', {'eid'});

    for p = 1:360
        brain_data_360p.(sprintf('SFC_P%03d', p)) = SCDFCC_Mean(:, p);
    end

    for p = 1:360
        brain_data_360p.(sprintf('DSFC_P%03d', p)) = SCDFCC_SD(:, p);
    end

    % ---------------------------------------------------------------------
    % Save
    % ---------------------------------------------------------------------

    tag = sprintf('TR%d', w);

    outputs.(tag).brain_data_7n = brain_data_7n;
    outputs.(tag).brain_data_360p = brain_data_360p;

    save( ...
        fullfile(analysis_table_dir, ...
        sprintf('brain_data_%dTR.mat', w)), ...
        'brain_data_7n', ...
        'brain_data_360p', ...
        'SFC_7n', ...
        'DSFC_7n', ...
        'common_eid', ...
        '-v7.3');

    writetable( ...
        brain_data_7n, ...
        fullfile(analysis_table_dir, ...
        sprintf('brain_data_7n_%dTR.csv', w)));

    writetable( ...
        brain_data_360p, ...
        fullfile(analysis_table_dir, ...
        sprintf('brain_data_360p_%dTR.csv', w)));

    fprintf('Prepared %d-TR analysis tables: N = %d\n', ...
        w, height(brain_data_7n));

end

end


function Xnet = network_mean(Xparcel, network_labels, n_networks)
% NETWORK_MEAN
% Average parcel-level values within each network.

Xnet = nan(size(Xparcel, 1), n_networks);

for n = 1:n_networks
    idx = network_labels == n;
    Xnet(:, n) = mean(Xparcel(:, idx), 2, 'omitnan');
end

end
