%% 03_prepare_sfc_analysis_tables.m
% Convert parcel-level SFC/DSFC arrays into analysis-ready MATLAB and CSV
% tables. The script processes every window length listed in config.m.

if ~exist('config.m', 'file')
    error('Copy config_example.m to config.m and configure local paths.');
end
run('config.m');
addpath(fullfile(PROJECT_ROOT, 'utilities'));

if ~exist(YEO7_GLASSER_FILE, 'file')
    error('Yeo-7 Glasser mapping not found: %s', YEO7_GLASSER_FILE);
end
load(YEO7_GLASSER_FILE, 'Yeo7MMP');
Yeo7MMP = Yeo7MMP(:);
if numel(Yeo7MMP) ~= N_CORTICAL_PARCELS
    error('Yeo7MMP must contain one label per Glasser parcel.');
end

network_names = {'Visual', 'Somatomotor', 'DorsalAttention', ...
    'VentralAttention', 'Limbic', 'Frontoparietal', 'Default'};
analysis_table_dir = fullfile(OUTPUT_DIR, 'analysis_tables');
if ~exist(analysis_table_dir, 'dir')
    mkdir(analysis_table_dir);
end

for window_length = WINDOW_LENGTHS
    input_file = fullfile(OUTPUT_DIR, sprintf( ...
        'ukb_glasser_sfc_dsfc_streamline_count_%dTR.mat', window_length));
    if ~exist(input_file, 'file')
        warning('Skipping missing input: %s', input_file);
        continue;
    end

    load(input_file, 'SCDFCC_Mean', 'SCDFCC_SD', 'common_eid');
    if size(SCDFCC_Mean, 2) ~= N_CORTICAL_PARCELS || ...
            size(SCDFCC_SD, 2) ~= N_CORTICAL_PARCELS
        error('Unexpected parcel count in %s.', input_file);
    end
    if size(SCDFCC_Mean, 1) ~= numel(common_eid) || ...
            size(SCDFCC_SD, 1) ~= numel(common_eid)
        error('Participant IDs and coupling arrays differ in %s.', input_file);
    end

    SFC_7n = nan(numel(common_eid), 7);
    DSFC_7n = nan(numel(common_eid), 7);
    for network_index = 1:7
        parcel_mask = Yeo7MMP == network_index;
        SFC_7n(:, network_index) = mean( ...
            SCDFCC_Mean(:, parcel_mask), 2, 'omitnan');
        DSFC_7n(:, network_index) = mean( ...
            SCDFCC_SD(:, parcel_mask), 2, 'omitnan');
    end

    brain_data_7n = table(common_eid(:), 'VariableNames', {'eid'});
    for network_index = 1:7
        brain_data_7n.(sprintf('SFC_%s', network_names{network_index})) = ...
            SFC_7n(:, network_index);
        brain_data_7n.(sprintf('DSFC_%s', network_names{network_index})) = ...
            DSFC_7n(:, network_index);
    end

    brain_data_360p = table(common_eid(:), 'VariableNames', {'eid'});
    for parcel_index = 1:N_CORTICAL_PARCELS
        brain_data_360p.(sprintf('SFC_P%03d', parcel_index)) = ...
            SCDFCC_Mean(:, parcel_index);
        brain_data_360p.(sprintf('DSFC_P%03d', parcel_index)) = ...
            SCDFCC_SD(:, parcel_index);
    end

    mat_output = fullfile(analysis_table_dir, sprintf( ...
        'brain_data_%dTR.mat', window_length));
    save(mat_output, 'brain_data_7n', 'brain_data_360p', ...
        'SFC_7n', 'DSFC_7n', 'common_eid', '-v7.3');
    writetable(brain_data_7n, fullfile(analysis_table_dir, sprintf( ...
        'brain_data_7n_%dTR.csv', window_length)));
    writetable(brain_data_360p, fullfile(analysis_table_dir, sprintf( ...
        'brain_data_360p_%dTR.csv', window_length)));
    fprintf('Prepared %d-TR tables for %d participants.\n', ...
        window_length, height(brain_data_7n));
end
