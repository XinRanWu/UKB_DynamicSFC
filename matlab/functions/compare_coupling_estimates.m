function [summary, network_correlations] = compare_coupling_estimates(reference, alternative, yeo_labels, label, spin_indices)
%COMPARE_COUPLING_ESTIMATES Compare matched coupling definitions.
if nargin < 5, spin_indices = []; end
[eid, i_ref, i_alt] = intersect(reference.eid,alternative.eid,'stable');
measures = {'SFC','DSFC'};
summary_rows = cell(2,8);
network_rows = cell(14,5); row = 0;
for m = 1:2
    x = reference.(measures{m})(i_ref,:); y = alternative.(measures{m})(i_alt,:);
    within = nan(numel(eid),1);
    for i = 1:numel(eid)
        within(i) = corr(x(i,:)',y(i,:)','Type','Spearman','Rows','pairwise');
    end
    if isempty(spin_indices)
        group_rho = corr(mean(x,1,'omitnan')',mean(y,1,'omitnan')', ...
            'Type','Spearman','Rows','pairwise');
        group_spin_p = NaN;
    else
        [group_rho,group_spin_p] = spin_test_correlation( ...
            mean(x,1,'omitnan')',mean(y,1,'omitnan')',spin_indices);
    end
    summary_rows(m,:) = {label,measures{m},numel(eid),mean(within,'omitnan'), ...
        median(within,'omitnan'),group_rho,group_spin_p,alternative.edge_rule};
    x7 = aggregate_yeo7(x,yeo_labels); y7 = aggregate_yeo7(y,yeo_labels);
    for network = 1:7
        row = row+1;
        [r,p] = corr(x7(:,network),y7(:,network),'Type','Pearson','Rows','pairwise');
        network_rows(row,:) = {label,measures{m},network,r,p};
    end
end
summary = cell2table(summary_rows,'VariableNames', ...
    {'comparison','measure','N','mean_within_spearman','median_within_spearman', ...
     'group_spearman','group_spin_p','edge_rule'});
network_correlations = cell2table(network_rows,'VariableNames', ...
    {'comparison','measure','network','pearson_r','p'});
end
