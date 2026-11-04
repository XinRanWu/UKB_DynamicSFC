%% Prepare Yeo-7 SFC and DSFC phenotypes for PLINK 2 GWAS
cfg = load_config();
load(cfg.analysis_coupling_table_file,'coupling_7n');
basic = read_analysis_table(cfg.basic_file,'basic/MRI covariates');
head_motion = basic(:,{'eid','HeadMotion_2_0'});
gwas = innerjoin_many(head_motion,coupling_7n);
gwas(gwas.HeadMotion_2_0>0.20,:) = [];
gwas.HeadMotion_2_0 = [];
gwas.FID = gwas.eid;
gwas.IID = gwas.eid;
gwas = movevars(gwas,{'FID','IID'},'Before','eid');
gwas.eid = [];

output_dir = fullfile(cfg.output_root,'gwas');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
writetable(gwas,fullfile(output_dir,'coupling_phenotypes.tsv'), ...
    'FileType','text','Delimiter','\t');
