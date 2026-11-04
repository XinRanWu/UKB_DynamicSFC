function [rho, p_spin] = spin_test_correlation(reference_map, comparison_map, spin_indices)
%SPIN_TEST_CORRELATION Two-sided spatial permutation test for map similarity.
% spin_indices is parcels x permutations; each column indexes a rotated map.
reference_map = reference_map(:); comparison_map = comparison_map(:);
if size(spin_indices,1) ~= numel(reference_map)
    error('spin_indices must have one row per cortical parcel.');
end
rho = corr(reference_map,comparison_map,'Type','Spearman','Rows','pairwise');
null_rho = nan(size(spin_indices,2),1);
for permutation = 1:size(spin_indices,2)
    index = spin_indices(:,permutation);
    valid = isfinite(index) & index>=1 & index<=numel(comparison_map);
    rotated = nan(size(comparison_map));
    rotated(valid) = comparison_map(index(valid));
    null_rho(permutation) = corr(reference_map,rotated,'Type','Spearman','Rows','pairwise');
end
p_spin = (1 + sum(abs(null_rho)>=abs(rho),'omitnan')) / ...
    (1 + sum(isfinite(null_rho)));
end
