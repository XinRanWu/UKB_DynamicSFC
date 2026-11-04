%% Build parcel- and network-level coupling tables
cfg = load_config();
input_file = fullfile(cfg.output_root,'coupling','primary_glasser360_30TR.mat');
if ~exist(input_file,'file'), error('Run 01_compute_primary_coupling.m first.'); end
load(input_file,'primary');
yeo = load_mat_variable(cfg.yeo7_glasser_file,cfg.yeo7_variable);
[coupling_360,coupling_7n] = coupling_to_tables(primary,yeo,cfg.network_names);

output_dir = fullfile(cfg.output_root,'tables');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
writetable(coupling_360,fullfile(output_dir,'coupling_360.csv'));
writetable(coupling_7n,fullfile(output_dir,'coupling_7n.csv'));
save(fullfile(output_dir,'coupling_tables.mat'),'coupling_360','coupling_7n','-v7.3');
fprintf('Analysis tables written for %d participants.\n',height(coupling_7n));
