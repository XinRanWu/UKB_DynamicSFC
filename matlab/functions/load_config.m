function cfg = load_config()
%LOAD_CONFIG Load the private MATLAB configuration from repository root.

this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(fileparts(this_file)));
config_file = fullfile(repo_root, 'config', 'config.m');
if ~exist(config_file, 'file')
    error('Missing config/config.m. Copy config/config.example.m first.');
end
run(config_file);
required = {'project_root','data_root','output_root','window_length', ...
    'window_step','n_parcels','network_names'};
for i = 1:numel(required)
    if ~isfield(cfg, required{i})
        error('Configuration field cfg.%s is required.', required{i});
    end
end
end
