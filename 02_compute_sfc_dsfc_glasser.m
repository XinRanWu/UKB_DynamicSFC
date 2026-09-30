%% 02_compute_sfc_dsfc_glasser.m
% Compute participant-level static SFC and dynamic SFC variability (DSFC).
% This is a script, not a function. Configure paths and parameters in config.m.

if ~exist('config.m', 'file')
    error('Copy config_example.m to config.m and configure local paths.');
end
has_window_override = exist('WINDOW_LENGTH', 'var') == 1;
has_step_override = exist('WINDOW_STEP', 'var') == 1;
has_worker_override = exist('N_WORKERS', 'var') == 1;
if has_window_override
    requested_window_length = WINDOW_LENGTH;
end
if has_step_override
    requested_window_step = WINDOW_STEP;
end
if has_worker_override
    requested_workers = N_WORKERS;
end
run('config.m');
if has_window_override
    WINDOW_LENGTH = requested_window_length;
end
if has_step_override
    WINDOW_STEP = requested_window_step;
end
if has_worker_override
    N_WORKERS = requested_workers;
end

required_files = {FMRI_FILE, SC_FILE};
for file_index = 1:numel(required_files)
    if ~exist(required_files{file_index}, 'file')
        error('Required input not found: %s', required_files{file_index});
    end
end

fprintf('Loading fMRI time series from %s\n', FMRI_FILE);
load(FMRI_FILE, 'fMRI_Glasser_TS_2_0', 'fMRI_Glasser_eID_2_0');
fprintf('Loading structural connectomes from %s\n', SC_FILE);
load(SC_FILE, 'connectome_streamline_count_10M_2_0', ...
    'dMRI_Glasser_S1_eID_2_0');

common_eid = intersect(fMRI_Glasser_eID_2_0, ...
    dMRI_Glasser_S1_eID_2_0, 'stable');
[found_fmri, fmri_order] = ismember(common_eid, fMRI_Glasser_eID_2_0);
[found_dmri, dmri_order] = ismember(common_eid, dMRI_Glasser_S1_eID_2_0);
if ~all(found_fmri) || ~all(found_dmri)
    error('Participant matching failed.');
end
fmri_data = fMRI_Glasser_TS_2_0(fmri_order);
sc_data = connectome_streamline_count_10M_2_0(dmri_order);
fprintf('Matched %d participants across modalities.\n', numel(common_eid));

parcel_pair_size = [];
if NORMALIZE_STREAMLINES_BY_PARCEL_SIZE
    if ~exist(PARCEL_SIZE_FILE, 'file')
        error(['Parcel-size normalization was requested, but the parcel-size ', ...
            'file is missing: %s'], PARCEL_SIZE_FILE);
    end
    parcel_size_data = load(PARCEL_SIZE_FILE);
    if ~isfield(parcel_size_data, 'roi_size') || ...
            numel(parcel_size_data.roi_size) ~= N_CORTICAL_PARCELS
        error('PARCEL_SIZE_FILE must contain a %d-element roi_size vector.', ...
            N_CORTICAL_PARCELS);
    end
    roi_size = parcel_size_data.roi_size(:);
    parcel_pair_size = roi_size * roi_size';
end

pool_object = gcp('nocreate');
if isempty(pool_object) && N_WORKERS > 0
    parpool('local', N_WORKERS);
end

n_participants = numel(common_eid);
SCDFCC_Mean = nan(n_participants, N_CORTICAL_PARCELS);
SCDFCC_SD = nan(n_participants, N_CORTICAL_PARCELS);
analysis_start = tic;

parfor participant_index = 1:n_participants
    structural_connectivity = sc_data{participant_index};
    structural_connectivity = structural_connectivity( ...
        1:N_CORTICAL_PARCELS, 1:N_CORTICAL_PARCELS);
    if NORMALIZE_STREAMLINES_BY_PARCEL_SIZE
        structural_connectivity = structural_connectivity ./ parcel_pair_size;
    end

    time_series = fmri_data{participant_index};
    if size(time_series, 1) == N_CORTICAL_PARCELS
        time_series = time_series;
    elseif size(time_series, 2) == N_CORTICAL_PARCELS
        time_series = time_series';
    else
        error('Participant %d fMRI data do not contain %d parcels.', ...
            participant_index, N_CORTICAL_PARCELS);
    end

    n_timepoints = size(time_series, 2);
    window_starts = 1:WINDOW_STEP:(n_timepoints - WINDOW_LENGTH + 1);
    if isempty(window_starts)
        error('WINDOW_LENGTH exceeds the available fMRI time series.');
    end

    window_sfc = nan(numel(window_starts), N_CORTICAL_PARCELS);
    for window_index = 1:numel(window_starts)
        first_timepoint = window_starts(window_index);
        timepoint_index = first_timepoint:(first_timepoint + WINDOW_LENGTH - 1);
        functional_connectivity = corr( ...
            time_series(:, timepoint_index)', 'Rows', 'pairwise');

        for parcel_index = 1:N_CORTICAL_PARCELS
            other_parcels = [1:(parcel_index - 1), ...
                (parcel_index + 1):N_CORTICAL_PARCELS];
            window_sfc(window_index, parcel_index) = corr( ...
                structural_connectivity(parcel_index, other_parcels)', ...
                functional_connectivity(parcel_index, other_parcels)', ...
                'Type', 'Spearman', 'Rows', 'pairwise');
        end
    end

    SCDFCC_Mean(participant_index, :) = mean(window_sfc, 1, 'omitnan');
    SCDFCC_SD(participant_index, :) = std(window_sfc, 0, 1, 'omitnan');
end

if EXCLUDE_INCOMPLETE_COUPLING
    complete_participant = all(isfinite(SCDFCC_Mean), 2) & ...
        all(isfinite(SCDFCC_SD), 2);
    fprintf('Excluding %d participants with incomplete coupling measures.\n', ...
        sum(~complete_participant));
    common_eid = common_eid(complete_participant);
    SCDFCC_Mean = SCDFCC_Mean(complete_participant, :);
    SCDFCC_SD = SCDFCC_SD(complete_participant, :);
end

if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end
output_file = fullfile(OUTPUT_DIR, sprintf( ...
    'ukb_glasser_sfc_dsfc_streamline_count_%dTR.mat', WINDOW_LENGTH));
save(output_file, 'common_eid', 'SCDFCC_Mean', 'SCDFCC_SD', ...
    'WINDOW_LENGTH', 'WINDOW_STEP', ...
    'NORMALIZE_STREAMLINES_BY_PARCEL_SIZE', ...
    'EXCLUDE_INCOMPLETE_COUPLING', '-v7.3');

fprintf('Saved %s\n', output_file);
fprintf('Elapsed time: %.2f seconds.\n', toc(analysis_start));
