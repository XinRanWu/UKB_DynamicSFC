function result = compute_sfc_dsfc(fmri_cells, fmri_eid, sc_cells, sc_eid, options)
%COMPUTE_SFC_DSFC Compute regional SFC and temporal variability of SFC.
%
% SFC is the mean across windows of the regional Spearman correlation
% between SC and FC profiles. DSFC is the standard deviation of the same
% window-wise SFC series. Self-connections are always excluded.
%
% options fields:
%   n_parcels, window_length, window_step, edge_rule, n_workers,
%   exclude_incomplete. edge_rule is all, positive_fc, nonzero_sc, or both.

arguments
    fmri_cells cell
    fmri_eid
    sc_cells cell
    sc_eid
    options.n_parcels (1,1) double {mustBeInteger,mustBePositive}
    options.window_length (1,1) double {mustBeInteger,mustBePositive}
    options.window_step (1,1) double {mustBeInteger,mustBePositive}
    options.edge_rule (1,:) char = 'all'
    options.n_workers (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.exclude_incomplete (1,1) logical = true
end

valid_rules = {'all','positive_fc','nonzero_sc','both'};
if ~ismember(options.edge_rule, valid_rules)
    error('edge_rule must be one of: %s.', strjoin(valid_rules, ', '));
end

fmri_eid = fmri_eid(:);
sc_eid = sc_eid(:);
[common_eid, fmri_index, sc_index] = intersect(fmri_eid, sc_eid, 'stable');
if isempty(common_eid)
    error('No participants are shared by the fMRI and SC inputs.');
end
fmri_cells = fmri_cells(fmri_index);
sc_cells = sc_cells(sc_index);

if options.n_workers > 0 && isempty(gcp('nocreate'))
    parpool('local', options.n_workers);
end

n = numel(common_eid);
sfc = nan(n, options.n_parcels);
dsfc = nan(n, options.n_parcels);

parfor participant = 1:n
    ts = orient_timeseries(fmri_cells{participant}, options.n_parcels);
    sc = sc_cells{participant};
    if size(sc,1) < options.n_parcels || size(sc,2) < options.n_parcels
        error('SC matrix for participant %d has fewer than %d parcels.', ...
            participant, options.n_parcels);
    end
    sc = double(sc(1:options.n_parcels, 1:options.n_parcels));

    starts = 1:options.window_step:(size(ts,2)-options.window_length+1);
    if isempty(starts)
        error('Window length exceeds the time series for participant %d.', participant);
    end
    window_sfc = nan(numel(starts), options.n_parcels);
    for w = 1:numel(starts)
        idx = starts(w):(starts(w)+options.window_length-1);
        fc = corr(ts(:,idx)', 'Type', 'Pearson', 'Rows', 'pairwise');
        for parcel = 1:options.n_parcels
            keep = true(1, options.n_parcels);
            keep(parcel) = false;
            switch options.edge_rule
                case 'positive_fc'
                    keep = keep & fc(parcel,:) > 0;
                case 'nonzero_sc'
                    keep = keep & sc(parcel,:) ~= 0;
                case 'both'
                    keep = keep & fc(parcel,:) > 0 & sc(parcel,:) ~= 0;
            end
            valid = keep & isfinite(fc(parcel,:)) & isfinite(sc(parcel,:));
            if nnz(valid) >= 3
                window_sfc(w,parcel) = corr(sc(parcel,valid)', ...
                    fc(parcel,valid)', 'Type', 'Spearman', 'Rows', 'complete');
            end
        end
    end
    sfc(participant,:) = mean(window_sfc, 1, 'omitnan');
    dsfc(participant,:) = std(window_sfc, 0, 1, 'omitnan');
end

if options.exclude_incomplete
    complete = all(isfinite(sfc),2) & all(isfinite(dsfc),2);
    common_eid = common_eid(complete);
    sfc = sfc(complete,:);
    dsfc = dsfc(complete,:);
end

result = struct();
result.eid = common_eid;
result.SFC = sfc;
result.DSFC = dsfc;
result.window_length = options.window_length;
result.window_step = options.window_step;
result.edge_rule = options.edge_rule;
result.n_parcels = options.n_parcels;
end

function ts = orient_timeseries(ts, n_parcels)
ts = double(ts);
if size(ts,1) == n_parcels
    return;
elseif size(ts,2) == n_parcels
    ts = ts';
else
    error('Time series has dimensions %d x %d; expected %d parcels.', ...
        size(ts,1), size(ts,2), n_parcels);
end
end
