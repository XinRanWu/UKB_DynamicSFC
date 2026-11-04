function [parcel_table,network_table] = bold_properties_to_tables(result,yeo_labels,network_names)
%BOLD_PROPERTIES_TO_TABLES Create named AC1 and BOLD-variability tables.
parcel_table = table(result.eid(:),'VariableNames',{'eid'});
for parcel = 1:size(result.AC1,2)
    parcel_table.(sprintf('AC1_P%03d',parcel)) = result.AC1(:,parcel);
    parcel_table.(sprintf('BV_P%03d',parcel)) = result.BOLD_variability(:,parcel);
end
ac7 = aggregate_yeo7(result.AC1,yeo_labels);
bv7 = aggregate_yeo7(result.BOLD_variability,yeo_labels);
network_table = table(result.eid(:),'VariableNames',{'eid'});
for network = 1:7
    network_table.(['AC1_' network_names{network}]) = ac7(:,network);
    network_table.(['BV_' network_names{network}]) = bv7(:,network);
end
end
