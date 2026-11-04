function network_values = aggregate_yeo7(parcel_values, yeo_labels)
%AGGREGATE_YEO7 Average parcel values within the seven Yeo networks.
yeo_labels = yeo_labels(:);
if size(parcel_values,2) ~= numel(yeo_labels)
    error('Parcel array and Yeo label vector have different lengths.');
end
if ~all(ismember(unique(yeo_labels), 1:7))
    error('Yeo labels must contain integers from 1 to 7.');
end
network_values = nan(size(parcel_values,1), 7);
for network = 1:7
    network_values(:,network) = mean(parcel_values(:,yeo_labels==network), ...
        2, 'omitnan');
end
end
