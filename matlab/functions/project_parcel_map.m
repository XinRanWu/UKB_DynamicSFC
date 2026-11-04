function vertex_map = project_parcel_map(parcel_map, vertex_labels)
%PROJECT_PARCEL_MAP Expand parcel values to a surface label vector.
parcel_map = parcel_map(:); vertex_labels = vertex_labels(:);
vertex_map = nan(size(vertex_labels));
valid = isfinite(vertex_labels) & vertex_labels>=1 & vertex_labels<=numel(parcel_map);
vertex_map(valid) = parcel_map(vertex_labels(valid));
end
