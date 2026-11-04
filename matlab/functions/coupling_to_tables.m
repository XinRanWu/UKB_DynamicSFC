function [parcel_table, network_table] = coupling_to_tables(result, yeo_labels, network_names)
%COUPLING_TO_TABLES Convert coupling arrays to named participant tables.
if numel(network_names) ~= 7
    error('Exactly seven network names are required.');
end
parcel_table = table(result.eid(:), 'VariableNames', {'eid'});
for parcel = 1:size(result.SFC,2)
    parcel_table.(sprintf('SFC_P%03d',parcel)) = result.SFC(:,parcel);
    parcel_table.(sprintf('DSFC_P%03d',parcel)) = result.DSFC(:,parcel);
end

sfc7 = aggregate_yeo7(result.SFC, yeo_labels);
dsfc7 = aggregate_yeo7(result.DSFC, yeo_labels);
network_table = table(result.eid(:), 'VariableNames', {'eid'});
for network = 1:7
    network_table.(['SFC_' network_names{network}]) = sfc7(:,network);
    network_table.(['DSFC_' network_names{network}]) = dsfc7(:,network);
end
end
