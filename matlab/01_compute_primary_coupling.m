%% Compute primary Glasser-360 SFC and DSFC
cfg = load_config();

fmri = load_mat_variable(cfg.fmri_file, cfg.fmri_data_variable);
fmri_eid = load_mat_variable(cfg.fmri_file, cfg.fmri_id_variable);
sc = load_mat_variable(cfg.sc_streamline_file, cfg.sc_streamline_variable);
sc_eid = load_mat_variable(cfg.sc_streamline_file, cfg.sc_id_variable);

options = struct('n_parcels',cfg.n_parcels, ...
    'window_length',cfg.window_length,'window_step',cfg.window_step, ...
    'edge_rule',cfg.edge_rule,'n_workers',cfg.n_workers, ...
    'exclude_incomplete',cfg.exclude_incomplete_coupling);
primary = compute_sfc_dsfc(fmri,fmri_eid,sc,sc_eid,options);

output_dir = fullfile(cfg.output_root,'coupling');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
output_file = fullfile(output_dir,'primary_glasser360_30TR.mat');
save(output_file,'primary','options','-v7.3');
fprintf('Primary coupling saved to %s\n',output_file);
