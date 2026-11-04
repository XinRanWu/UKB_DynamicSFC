function result = compute_bold_properties(fmri_cells, fmri_eid, n_parcels)
%COMPUTE_BOLD_PROPERTIES Regional lag-1 autocorrelation and BOLD variability.
n = numel(fmri_eid); ac1 = nan(n,n_parcels); variability = nan(n,n_parcels);
parfor participant = 1:n
    ts = double(fmri_cells{participant});
    if size(ts,1) ~= n_parcels && size(ts,2) == n_parcels, ts = ts'; end
    if size(ts,1) ~= n_parcels, error('Unexpected fMRI matrix orientation.'); end
    variability(participant,:) = std(ts,0,2,'omitnan')';
    for parcel = 1:n_parcels
        ac1(participant,parcel) = corr(ts(parcel,1:end-1)',ts(parcel,2:end)', ...
            'Type','Pearson','Rows','pairwise');
    end
end
result = struct('eid',fmri_eid(:),'AC1',ac1,'BOLD_variability',variability);
end
