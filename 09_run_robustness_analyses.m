%% 09_run_robustness_analyses.m
% Quantify robustness across window lengths, structural-connectivity measures,
% and parcellations. Results are written as inspectable CSV tables.

if ~exist('config.m', 'file')
    error('Copy config_example.m to config.m and configure local paths.');
end
run('config.m');

primary_file = fullfile(OUTPUT_DIR, sprintf( ...
    'ukb_glasser_sfc_dsfc_streamline_count_%dTR.mat', WINDOW_LENGTH));
if ~exist(primary_file, 'file')
    error('Primary coupling file not found: %s', primary_file);
end
primary = load(primary_file, 'common_eid', 'SCDFCC_Mean', 'SCDFCC_SD');

summary_rows = {};
comparison_index = 0;
for alternate_window = setdiff(WINDOW_LENGTHS, WINDOW_LENGTH)
    alternate_file = fullfile(OUTPUT_DIR, sprintf( ...
        'ukb_glasser_sfc_dsfc_streamline_count_%dTR.mat', alternate_window));
    if ~exist(alternate_file, 'file')
        warning('Skipping missing window-length result: %s', alternate_file);
        continue;
    end
    alternate = load(alternate_file, ...
        'common_eid', 'SCDFCC_Mean', 'SCDFCC_SD');
    [matched_eid, primary_index, alternate_index] = intersect( ...
        primary.common_eid, alternate.common_eid, 'stable');

    measures = {'SFC', 'DSFC'};
    primary_arrays = {primary.SCDFCC_Mean, primary.SCDFCC_SD};
    alternate_arrays = {alternate.SCDFCC_Mean, alternate.SCDFCC_SD};
    for measure_index = 1:2
        primary_values = primary_arrays{measure_index}(primary_index, :);
        alternate_values = alternate_arrays{measure_index}(alternate_index, :);
        individual_r = nan(numel(matched_eid), 1);
        for participant_index = 1:numel(matched_eid)
            individual_r(participant_index) = corr( ...
                primary_values(participant_index, :)', ...
                alternate_values(participant_index, :)', ...
                'Type', 'Spearman', 'Rows', 'pairwise');
        end
        group_r = corr(mean(primary_values, 1, 'omitnan')', ...
            mean(alternate_values, 1, 'omitnan')', ...
            'Type', 'Spearman', 'Rows', 'pairwise');
        comparison_index = comparison_index + 1;
        summary_rows(comparison_index, :) = { ...
            sprintf('%dTR_vs_%dTR', WINDOW_LENGTH, alternate_window), ...
            measures{measure_index}, numel(matched_eid), ...
            mean(individual_r, 'omitnan'), median(individual_r, 'omitnan'), ...
            group_r};
    end
end

if ~isempty(FA_SFC_FILE)
    if ~exist(FA_SFC_FILE, 'file')
        warning('Skipping missing FA comparison file: %s', FA_SFC_FILE);
    else
        alternate = load(FA_SFC_FILE, ...
            'common_eid', 'SCDFCC_Mean', 'SCDFCC_SD');
        [matched_eid, primary_index, alternate_index] = intersect( ...
            primary.common_eid, alternate.common_eid, 'stable');
        measures = {'SFC', 'DSFC'};
        primary_arrays = {primary.SCDFCC_Mean, primary.SCDFCC_SD};
        alternate_arrays = {alternate.SCDFCC_Mean, alternate.SCDFCC_SD};
        for measure_index = 1:2
            primary_values = primary_arrays{measure_index}(primary_index, :);
            alternate_values = alternate_arrays{measure_index}(alternate_index, :);
            individual_r = nan(numel(matched_eid), 1);
            for participant_index = 1:numel(matched_eid)
                individual_r(participant_index) = corr( ...
                    primary_values(participant_index, :)', ...
                    alternate_values(participant_index, :)', ...
                    'Type', 'Spearman', 'Rows', 'pairwise');
            end
            group_r = corr(mean(primary_values, 1, 'omitnan')', ...
                mean(alternate_values, 1, 'omitnan')', ...
                'Type', 'Spearman', 'Rows', 'pairwise');
            comparison_index = comparison_index + 1;
            summary_rows(comparison_index, :) = { ...
                'streamline_count_vs_FA', measures{measure_index}, ...
                numel(matched_eid), mean(individual_r, 'omitnan'), ...
                median(individual_r, 'omitnan'), group_r};
        end
    end
end

if isempty(summary_rows)
    warning('No robustness comparisons were available.');
else
    robustness_summary = cell2table(summary_rows, 'VariableNames', ...
        {'comparison', 'measure', 'n', 'mean_individual_r', ...
         'median_individual_r', 'group_spatial_r'});
    robustness_dir = fullfile(OUTPUT_DIR, 'robustness_results');
    if ~exist(robustness_dir, 'dir')
        mkdir(robustness_dir);
    end
    output_file = fullfile(robustness_dir, 'robustness_summary.csv');
    writetable(robustness_summary, output_file);
    fprintf('Robustness summary saved to %s\n', output_file);
end

if ~isempty(SCHAEFER_SFC_FILE)
    warning(['Schaefer-to-Glasser spatial comparison requires a documented ', ...
        'surface projection. The original exploratory code remains in Git ', ...
        'history until those atlas resources can be redistributed or cited ', ...
        'with exact versions.']);
end
