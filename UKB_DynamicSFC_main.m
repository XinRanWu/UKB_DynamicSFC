load('/public/home/zhangjie/DataAnalysis/wxr_toolbox/toolbox_matlab/cbrewer2-master/cbrewer2/colorbrewer.mat', 'colorbrewer');
UKB_Outcomes_2_0 = readtable('/public/home/zhangjie/UKB_Outcomes_2_0.csv');
UKB_Lifestyle_0_0 = readtable('/public/home/zhangjie/UKB_Lifestyle_MRI_0_0.csv');
Yeo7Names = {'Visual';'Somatomotor';'DorsalAttention';'VentralAttention';'Limbic';'Frontoparietal';'Default'};

UKB_Lifestyle_0_0.x1160_L = double(UKB_Lifestyle_0_0.x1160 >= 9);
UKB_Lifestyle_0_0.x1160_S = double(UKB_Lifestyle_0_0.x1160 < 6);
UKB_Lifestyle_0_0.x1160_S(isnan(UKB_Lifestyle_0_0.x1160)) = NaN;
UKB_Lifestyle_0_0.x1160_L(isnan(UKB_Lifestyle_0_0.x1160)) = NaN;
UKB_Lifestyle_0_0 = removevars(UKB_Lifestyle_0_0, 'x1160');
x709 = UKB_Lifestyle_0_0.x709;
UKB_Lifestyle_0_0.x709 = double(UKB_Lifestyle_0_0.x709 == 1);
UKB_Lifestyle_0_0.x709(isnan(x709)) = NaN;
UKB_Lifestyle_0_0 = removevars(UKB_Lifestyle_0_0, {'x1289','x1299','x1309','x1319','x1329','x1339','x1349','x1369','x1379','x1389','x1438','x1448','x1458'});
UKB_Lifestyle_0_0 = removevars(UKB_Lifestyle_0_0, {'x1568','x1578','x1588','x1598','x1608'});
UKB_Lifestyle_0_0 = innerjoin(UKB_Lifestyle_0_0,UKB_Diet_0_0);
UKB_Outcomes_2_0 = removevars(UKB_Outcomes_2_0, 'x20012');

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/atlas/HCPMMP_atlas_info.mat', 'Yeo7MMP')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/ukb_glasser_scfcc_2_0_streamline_count_30TR.mat', 'SCDFCC_Mean', 'SCDFCC_SD', 'common_eid')
SFC_7n = netMean1D(SCDFCC_Mean,Yeo7MMP);
DSFC_7n = netMean1D(SCDFCC_SD,Yeo7MMP);

brain_data_7n = [array2table(common_eid),array2table(SFC_7n),array2table(DSFC_7n)];
brain_data_7n.Properties.VariableNames{1} = 'eid';
brain_data_360p = [array2table(common_eid),array2table(SCDFCC_Mean),array2table(SCDFCC_SD)];
brain_data_360p.Properties.VariableNames{1} = 'eid';

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/ukb_glasser_scfcc_2_0_streamline_count_40TR.mat', 'SCDFCC_Mean', 'SCDFCC_SD', 'common_eid')
SFC_7n_40TR = netMean1D(SCDFCC_Mean,Yeo7MMP);
DSFC_7n_40TR = netMean1D(SCDFCC_SD,Yeo7MMP);

brain_data_7n_40TR = [array2table(common_eid),array2table(SFC_7n_40TR),array2table(DSFC_7n_40TR)];
brain_data_7n_40TR.Properties.VariableNames{1} = 'eid';
brain_data_360p_40TR = [array2table(common_eid),array2table(SCDFCC_Mean),array2table(SCDFCC_SD)];
brain_data_360p_40TR.Properties.VariableNames{1} = 'eid';

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/ukb_glasser_scfcc_2_0_streamline_count_50TR.mat', 'SCDFCC_Mean', 'SCDFCC_SD', 'common_eid')
SFC_7n_50TR = netMean1D(SCDFCC_Mean,Yeo7MMP);
DSFC_7n_50TR = netMean1D(SCDFCC_SD,Yeo7MMP);

brain_data_7n_50TR = [array2table(common_eid),array2table(SFC_7n_40TR),array2table(DSFC_7n_50TR)];
brain_data_7n_50TR.Properties.VariableNames{1} = 'eid';
brain_data_360p_50TR = [array2table(common_eid),array2table(SCDFCC_Mean),array2table(SCDFCC_SD)];
brain_data_360p_50TR.Properties.VariableNames{1} = 'eid';

% 7 networks
UKB_BrainXY = innerjoin(innerjoin(UKB_Basic,UKB_BrainMRICov),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; 

formula = ['AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0'};
results_age_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
UKB_BrainXY.AgeAttend_2_0_Sq = zscore(UKB_BrainXY.AgeAttend_2_0).^2;
formula = ['AgeAttend_2_0+Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0_Sq'};
results_age2_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+SNR_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+t1_vol_scaling_2_0+(1|Centre_0_0)'];
XList = {'Sex_0_0'};
results_sex_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+HeadMotion_2_0+SNR_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+t1_vol_scaling_2_0+(1|Centre_0_0)'];
XList = {'BMI_2_0'};
results_bmi_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+SNR_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+(1|Centre_0_0)'];
XList = {'t1_vol_scaling_2_0'};
results_tiv_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
results_demo_7n = [results_age_7n;results_age2_7n;results_sex_7n;results_bmi_7n;results_tiv_7n];

X1List = {'x20016','x137'};
X2List = {'AgeAttend_2_0'};
formula = ['Sex_0_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+SNR_2_0+',...
    'Vol_WB_TIV_2_0+Race_1+Race_2+Race_3+Race_4+(1|Centre_0_0)'];
results_interact_7n = LMM_withFormula_Interact(UKB_BrainXY,formula,X1List,X2List,YList);

UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_Basic,UKB_BrainMRICov),UKB_Outcomes_2_0),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; 
XList = UKB_Outcomes_2_0.Properties.VariableNames(2:end)';
formula = ['AgeAttend_2_0 + AgeAttend_2_0^2 + Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + SNR_2_0 + ',...
    'Race_1 + Race_2 + Race_3 + Race_4 + t1_vol_scaling_2_0 + (1|Centre_0_0)'];
results_pheno_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_Basic,UKB_BrainMRICov),UKB_Lifestyle_0_0(:,[1,95:end])),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; 
XList = UKB_Lifestyle_0_0.Properties.VariableNames(95:end)';
formula = ['AgeAttend_2_0 + AgeAttend_2_0^2 + Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + SNR_2_0 + ',...
    'Race_1 + Race_2 + Race_3 + Race_4 + t1_vol_scaling_2_0 + (1|Centre_0_0)'];
results_lifestyle_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);


UKB_PRS = readtable('/home1/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_PRS.tsv', 'FileType', 'text');
UKB_PRS.p26290 = str2double(UKB_PRS.p26290);

APOE4_genetype = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/gene/APOE4_genetype.txt');
APOE4_genetype = removevars(APOE4_genetype, 'FID');
APOE4_genetype = removevars(APOE4_genetype, {'PAT','MAT','SEX','PHENOTYPE'});

UKB_GeneCov = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/gene/UKB_GeneCov.txt','Delimiter' ,',','TreatAsEmpty','NA');
UKB_GeneCov.Properties.VariableNames{1} = 'eid';

% UKB_PRS = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_PRS.csv');
% UKB_PRS.venous_thromboembolic_disease_VTE_1 = str2double(UKB_PRS.venous_thromboembolic_disease_VTE_1);

UKB_PRS = readtable('/home1/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_PRS.tsv', 'FileType', 'text');
UKB_PRS.p26290 = str2double(UKB_PRS.p26290);

UKB_BrainXY = innerjoin(innerjoin(innerjoin(innerjoin(UKB_Basic,UKB_BrainMRICov),UKB_GeneCov),UKB_PRS),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = [brain_data_7n.Properties.VariableNames(2:end)']; 
XList = UKB_PRS.Properties.VariableNames(2:end)';
formula = ['AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+SNR_2_0+',...
    't1_vol_scaling_2_0+x22000_0_0+x22009_0_1+x22009_0_2+x22009_0_3+',...
    'x22009_0_4+x22009_0_5+x22009_0_6+x22009_0_7+x22009_0_8+x22009_0_9+x22009_0_10+(1|Centre_0_0)'];
results_PRS_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

UKB_BrainXY = innerjoin(innerjoin(innerjoin(innerjoin(UKB_Basic,UKB_BrainMRICov),UKB_GeneCov),APOE4_genetype),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';

XList = {'APOE4_dosage'};
formula = ['AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+SNR_2_0+',...
    't1_vol_scaling_2_0+x22000_0_0+x22009_0_1+x22009_0_2+x22009_0_3+',...
    'x22009_0_4+x22009_0_5+x22009_0_6+x22009_0_7+x22009_0_8+x22009_0_9+x22009_0_10+(1|Centre_0_0)'];

results_APOE4_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

AD_genetype = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/gene/AD_genetype.txt','Delimiter' ,',','TreatAsEmpty','NA');
AD_genetype = removevars(AD_genetype, 'FID');
AD_genetype = removevars(AD_genetype, {'PAT','MAT','SEX','PHENOTYPE'});

UKB_BrainXY = innerjoin(innerjoin(innerjoin(innerjoin(UKB_Basic,UKB_BrainMRICov),UKB_GeneCov),AD_genetype),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];

YList = brain_data_7n.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';
XList = AD_genetype.Properties.VariableNames(2:end)'; 
formula = ['AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+SNR_2_0+',...
    't1_vol_scaling_2_0+x22000_0_0+x22009_0_1+x22009_0_2+x22009_0_3+',...
    'x22009_0_4+x22009_0_5+x22009_0_6+x22009_0_7+x22009_0_8+x22009_0_9+x22009_0_10+(1|Centre_0_0)'];

results_AD_genes_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

X1List = {'APOE4_dosage'};
X2List = {'AgeAttend_2_0'};
formula = ['Sex_0_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+SNR_2_0+',...
    't1_vol_scaling_2_0+x22000_0_0+x22009_0_1+x22009_0_2+x22009_0_3+',...
    'x22009_0_4+x22009_0_5+x22009_0_6+x22009_0_7+x22009_0_8+x22009_0_9+x22009_0_10+(1|Centre_0_0)'];
results_GeneInteract_7n = LMM_withFormula_Interact(UKB_BrainXY,formula,X1List,X2List,YList);



X1List = {'x20016','x137'};
X2List = {'AgeAttend_2_0'};
formula = ['Sex_0_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+SNR_2_0+',...
    't1_vol_scaling_2_0+Race_1+Race_2+Race_3+Race_4+(1|Centre_0_0)'];
results_interact_360p = LMM_withFormula_Interact(UKB_BrainXY,formula,X1List,X2List,YList);


UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_Basic,UKB_BrainMRICov),UKB_Outcomes_2_0),brain_data_360p);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_360p.Properties.VariableNames(2:end)'; 
XList = {'x135','x137','x1970','x2010','x2178','x2188','x4282','x4526','x4548','x6373','x20016','x20197'};
formula = ['AgeAttend_2_0 + AgeAttend_2_0^2 + Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + SNR_2_0 + ',...
    'Race_1 + Race_2 + Race_3 + Race_4 + t1_vol_scaling_2_0 + (1|Centre_0_0)'];
results_pheno_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);



UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_Basic,UKB_BrainMRICov),UKB_Lifestyle_0_0(:,[1,95:end])),brain_data_360p);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_360p.Properties.VariableNames(2:end)'; 
XList = {'x738','x20022','x20116','x20495','x20497','x22189','x1160_L'}
formula = ['AgeAttend_2_0 + AgeAttend_2_0^2 + Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + SNR_2_0 + ',...
    'Race_1 + Race_2 + Race_3 + Race_4 + t1_vol_scaling_2_0 + (1|Centre_0_0)'];
results_lifestyle_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);


UKB_BrainXY = innerjoin(innerjoin(innerjoin(innerjoin(UKB_Basic,UKB_BrainMRICov),UKB_GeneCov),UKB_PRS),brain_data_360p);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = [brain_data_360p.Properties.VariableNames(2:end)']; 
XList = {'p26206','p26275'};
formula = ['AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+SNR_2_0+',...
    't1_vol_scaling_2_0+x22000_0_0+x22009_0_1+x22009_0_2+x22009_0_3+',...
    'x22009_0_4+x22009_0_5+x22009_0_6+x22009_0_7+x22009_0_8+x22009_0_9+x22009_0_10+(1|Centre_0_0)'];
results_PRS_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

UKB_BrainXY = innerjoin(innerjoin(innerjoin(innerjoin(UKB_Basic,UKB_BrainMRICov),UKB_GeneCov),APOE4_genetype),brain_data_360p);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_360p.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';

XList = {'APOE4_dosage'};
formula = ['AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+SNR_2_0+',...
    't1_vol_scaling_2_0+x22000_0_0+x22009_0_1+x22009_0_2+x22009_0_3+',...
    'x22009_0_4+x22009_0_5+x22009_0_6+x22009_0_7+x22009_0_8+x22009_0_9+x22009_0_10+(1|Centre_0_0)'];

results_APOE4_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);



X1List = {'APOE4_dosage'};
X2List = {'AgeAttend_2_0'};

formula = ['Sex_0_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+SNR_2_0+',...
    't1_vol_scaling_2_0+x22000_0_0+x22009_0_1+x22009_0_2+x22009_0_3+',...
    'x22009_0_4+x22009_0_5+x22009_0_6+x22009_0_7+x22009_0_8+x22009_0_9+x22009_0_10+(1|Centre_0_0)'];
results_GeneInteract_360p = LMM_withFormula_Interact(UKB_BrainXY,formula,X1List,X2List,YList);

%% GWAS
% load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/ukb_glasser_scfcc_2_0_mean_FA.mat', 'common_eid')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/ukb_glasser_scfcc_2_0_streamline_count.mat','common_eid', 'SCDFCC_Mean', 'SCDFCC_SD')
SFC_7n = netMean1D(SCDFCC_Mean,Yeo7MMP);
DSFC_7n = netMean1D(SCDFCC_SD,Yeo7MMP);
SFC_7n_Total = nanmean(SCDFCC_Mean,2);
DSFC_7n_Total = nanmean(SCDFCC_SD,2);
brain_data_7n = [array2table(common_eid),array2table([SFC_7n,SFC_7n_Total,DSFC_7n,DSFC_7n_Total])];
brain_data_7n.Properties.VariableNames{1} = 'eid';

brain_data_7n_tmp = innerjoin(UKB_Lifestyle_0_0(:,[1,77]),brain_data_7n);
brain_data_7n_tmp(brain_data_7n_tmp.HeadMotion_2_0>0.2,:) = [];
brain_data_7n_tmp = removevars(brain_data_7n_tmp, 'HeadMotion_2_0');

nx = {'SFC_VN';'SFC_SMN';'SFC_DAN';'SFC_VAN';'SFC_LN';'SFC_FPN';'SFC_DMN';'SFC_Global';...
    'DSFC_VN';'DSFC_SMN';'DSFC_DAN';'DSFC_VAN';'DSFC_LN';'DSFC_FPN';'DSFC_DMN';'DSFC_Global'};
brain_data_7n_tmp.Properties.VariableNames(2:end) = nx;


Xpheno = [array2table(brain_data_7n_tmp.eid),brain_data_7n_tmp];
Xpheno.Properties.VariableNames{2} = 'IID';Xpheno.Properties.VariableNames{1} = 'FID';
writetable(Xpheno,'/public/home/zhangjie/ZJLab/UKBiobank_Project/data/gene/UKB_dsfc_pheno_30TRfn.txt','Delimiter' ,' ')

cd('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/gene/')
phenofile = '/public/home/zhangjie/ZJLab/UKBiobank_Project/data/gene/UKB_dsfc_pheno_30TRfn.txt';

covarfile = '/public/home/zhangjie/ZJLab/UKBiobank_Project/data/gene/UKB_GeneCov_MRI.txt';
covarname = ['Sex_0_0,AgeAttend_2_0,BMI_2_0,HeadMotion_2_0,Vol_WB_TIV_2_0,Centre_2_0,',...
    'PC1-PC20'];


PList = Xpheno.Properties.VariableNames(3:end)';
parfor chr = 1:22
    bfile_qc = ['/public/home/ISTBI_data/UKB/gene/v3/QCed/ukb_imp_chr',num2str(chr),'_v3'];
    outfile = ['/public/home/zhangjie/ZJLab/UKBiobank_Project/project/',...
        'Structure_Function_Coupling_Aging/GWAS_results/gwas_chr',num2str(chr),'_7n'];
    plink_command = ...
        ['/public/home/zhangjie/DataAnalysis/wxr_toolbox/toobox_other/plink2_linux_x86_64_20191030/plink2 ',...
        '--bfile ',bfile_qc,' ',...
        '--linear  hide-covar --debug ',...
        '--pheno ',phenofile,' ',...
        '--pheno-name ',strjoin(PList, ','),' '...
        '--covar ',covarfile,' ',...
        '--covar-name ',covarname,' ',...
        '--covar-variance-standardize ',...
        '--out ',outfile];
    unix(plink_command);
    disp(['GWAS of in chr',num2str(chr)])
end
cd('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/GWAS_results/')

parfor i = 1:size(PList,1)
    
    unix(['head -n 1 gwas_chr1_7n.',PList{i},'.glm.linear > gwas.',PList{i},'.glm.linear']);
    unix(['tail -n +2 gwas_chr1_7n.',PList{i},'.glm.linear >> gwas.',PList{i},'.glm.linear']);
    
    for chr = 2:22
        unix(['tail -n +2 gwas_chr',num2str(chr),'_7n.',PList{i},'.glm.linear >> gwas.',PList{i},'.glm.linear']);
    end
    
    unix(sprintf(['sed -i "1s/.*/CHR\tPOS\tID\tREF\tALT\tPROVISIONAL_REF\tA1\tOMITTED\tA1_FREQ\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP\tERRCODE/" gwas.%s.glm.linear'], PList{i}));

    unix(['gzip gwas.',PList{i},'.glm.linear'])
end

