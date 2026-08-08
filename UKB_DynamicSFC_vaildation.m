UKB_Basic = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_Basic.csv','TreatAsEmpty','NA');
UKB_BrainMRICov = readtable('/home1/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_BrainMRICov.tsv', 'FileType', 'text');
UKB_BrainMRICov.t1_R_HBF_3_0 = str2double(UKB_BrainMRICov.t1_R_HBF_3_0);
UKB_BrainMRICov(sum(isnan(UKB_BrainMRICov{:,2:end}),2)==26,:) = [];

% UKB_BrainMRICov.Properties.VariableNames{22} = 't1_vol_scaling_2_0';
UKB_BrainMRICov.Properties.VariableNames{23} = 'Vol_WB_TIV_3_0';
UKB_BrainMRICov.Properties.VariableNames{24} = 'HBF_L_HBF_2_0';
UKB_BrainMRICov.Properties.VariableNames{26} = 'HBF_R_HBF_2_0';
UKB_BrainMRICov.Properties.VariableNames{25} = 'HBF_L_HBF_3_0';
UKB_BrainMRICov.Properties.VariableNames{27} = 'HBF_R_HBF_3_0';
UKB_BrainMRICov.Properties.VariableNames{14} = 'HeadMotion_2_0';
UKB_BrainMRICov.Properties.VariableNames{15} = 'HeadMotion_3_0';
UKB_BrainMRICov.Properties.VariableNames{8} = 'SNR_2_0';
UKB_BrainMRICov.Properties.VariableNames{9} = 'SNR_3_0';
UKB_BrainMRICov.Euler_Num_L_2_0 = 2 - 2 .* UKB_BrainMRICov.HBF_L_HBF_2_0;
UKB_BrainMRICov.Euler_Num_R_2_0 = 2 - 2 .* UKB_BrainMRICov.HBF_R_HBF_2_0;
UKB_BrainMRICov.Euler_Num_WB_2_0 = UKB_BrainMRICov.Euler_Num_R_2_0 + UKB_BrainMRICov.Euler_Num_L_2_0;


load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/atlas/HCPMMP_atlas_info.mat', 'Yeo7MMP')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/ukb_glasser_scfcc_2_0_streamline_count_30TR.mat', 'SCDFCC_Mean', 'SCDFCC_SD', 'common_eid')
SFC_7n = netMean1D(SCDFCC_Mean,Yeo7MMP);
DSFC_7n = netMean1D(SCDFCC_SD,Yeo7MMP);

load('/public/home/zhangjie/ZJLab/ZJLab_Toolset/rotate_parcellation-master/perm_centroid_info_HCPMMP1.mat')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/atlas/HCPMMP_atlas_info.mat', 'Yeo7MMP')

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

parfor i = 1:size(brain_data_7n,1)
    [r_SFC(i,1),p_SFC(i,1)] = corr(brain_data_360p{i,2:361}',brain_data_360p_40TR{i,2:361}');
    [r_DSFC(i,1),p_DSFC(i,1)] = corr(brain_data_360p{i,362:721}',brain_data_360p_40TR{i,362:721}');
    disp(num2str(i))
end

[r_7n, p_7n] = corr(brain_data_7n{:,2:15}, brain_data_7n_40TR{:,2:15},'Rows','pairwise');

id = ~isnan(r_SFC);
nanmean(r_SFC(id))
mean(p_SFC(id)<0.05)

id = ~isnan(r_DSFC);
nanmean(r_DSFC(id))
mean(p_DSFC(id)<0.05)

Yeo7Colormap = [0.471,0.0710,0.522;0.275,0.510,0.706;0,0.463,0.0550;...
    0.769,0.224,0.976;0.863,0.973,0.639;0.902,0.576,0.129;0.804,0.239,0.306];
Yeo7Names2 = {'VN','SMN','DAN','VAN','LN','FPN','DMN'};

averaged_DSFC = nanmean(brain_data_360p{:,362:721},1);
averaged_SFC = nanmean(brain_data_360p{:,2:361},1);
averaged_DSFC_40TR = nanmean(brain_data_360p_40TR{:,362:721},1);
averaged_SFC_40TR = nanmean(brain_data_360p_40TR{:,2:361},1);

averaged_DSFC_50TR = nanmean(brain_data_360p_50TR{:,362:721},1);
averaged_SFC_50TR = nanmean(brain_data_360p_50TR{:,2:361},1);

figure;

% Subplot 1: Scatter plot with colormap and black fit line
subplot(2,3,1);
scatter(averaged_SFC, averaged_SFC_40TR, 10, Yeo7MMP, 'filled');
colormap(Yeo7Colormap);
title('SFC');
hold on;
p1 = polyfit(averaged_SFC, averaged_SFC_40TR, 1); % Linear fit
f1 = polyval(p1, averaged_SFC);
plot(averaged_SFC, f1, 'k-', 'LineWidth', 1); % Black fit line
hold off;

% Subplot 2: Histogram with black mean vertical line
subplot(2,3,2);
histogram(r_SFC, 50);
title('SFC');
hold on;
mean_val_SFC = nanmean(r_SFC);
xline(mean_val_SFC, 'k', 'LineWidth', 1); % Black vertical line at the mean
hold off;

% Subplot 4: Scatter plot with colormap and black fit line
subplot(2,3,4);
scatter(averaged_DSFC, averaged_DSFC_40TR, 10, Yeo7MMP, 'filled');
colormap(Yeo7Colormap);
hold on;
p2 = polyfit(averaged_DSFC, averaged_DSFC_40TR, 1); % Linear fit
f2 = polyval(p2, averaged_DSFC);
plot(averaged_DSFC, f2, 'k-', 'LineWidth', 1); % Black fit line
hold off;
title('DSFC');

% Subplot 5: Histogram with black mean vertical line
subplot(2,3,5);
histogram(r_DSFC, 50);
title('DSFC');
hold on;
mean_val_DSFC = nanmean(r_DSFC);
xline(mean_val_DSFC, 'k', 'LineWidth', 1); % Black vertical line at the mean
hold off;

% Subplot 3: Heatmap with custom colormap (using heatmap function)
subplot(2,3,3);
[corr_matrix, p_matrix] = corr(brain_data_7n{:,2:8}, brain_data_7n_40TR{:,2:8}, 'Rows', 'pairwise');
ICC = [diag(corr_matrix),diag(p_matrix)];
h3 = heatmap(Yeo7Names2, Yeo7Names2, corr_matrix, 'ColorbarVisible', 'on', 'Colormap', hot(100)); % Using heatmap function
h3.XLabel = 'SFC';
h3.YLabel = 'SFC';
h3.Title = 'Correlation Heatmap';

% Subplot 6: Heatmap with custom colormap (using heatmap function)
subplot(2,3,6);
[corr_matrix_2, p_matrix_2] = corr(brain_data_7n{:,9:15}, brain_data_7n_40TR{:,9:15}, 'Rows', 'pairwise');
ICC = [diag(corr_matrix_2),diag(p_matrix_2)];
h6 = heatmap(Yeo7Names2, Yeo7Names2, corr_matrix_2, 'ColorbarVisible', 'on', 'Colormap', hot(100)); % Using heatmap function
h6.XLabel = 'DSFC';
h6.YLabel = 'DSFC';
h6.Title = 'Correlation Heatmap';


load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/ukb_Scheafer7n200p_scfcc_2_0_streamline_count.mat', 'SCDFCC_Mean', 'SCDFCC_SD', 'common_eid')

Scheafer200p_7Network = readtable('/home1/zhangjie/ZJLab/UKBiobank_Project/data/imaging/atlas/Scheafer200p_7Network.csv');

SFC_7n_Scheafer = netMean1D(SCDFCC_Mean,Scheafer200p_7Network.NetworkID);
DSFC_7n_Scheafer = netMean1D(SCDFCC_SD,Scheafer200p_7Network.NetworkID);

brain_data_7n_Scheafer = [array2table(common_eid),array2table(SFC_7n_Scheafer),array2table(DSFC_7n_Scheafer)];
brain_data_7n_Scheafer.Properties.VariableNames{1} = 'eid';
brain_data_200p_Scheafer = [array2table(common_eid),array2table(SCDFCC_Mean),array2table(SCDFCC_SD)];
brain_data_200p_Scheafer.Properties.VariableNames{1} = 'eid';

[s2f,f2s] = idFinderNum(brain_data_360p.eid,brain_data_200p_Scheafer.eid);

[r_7n, p_7n] = corr(brain_data_7n{s2f,2:15}, brain_data_7n_Scheafer{f2s,2:15},'Rows','pairwise');


Yeo7Colormap = [0.471,0.0710,0.522;0.275,0.510,0.706;0,0.463,0.0550;...
    0.769,0.224,0.976;0.863,0.973,0.639;0.902,0.576,0.129;0.804,0.239,0.306];
Yeo7Names2 = {'VN','SMN','DAN','VAN','LN','FPN','DMN'};

averaged_DSFC = nanmean(brain_data_360p{:,362:721},1);
averaged_SFC = nanmean(brain_data_360p{:,2:361},1);
averaged_DSFC_Scheafer = nanmean(brain_data_200p_Scheafer{:,202:401},1);
averaged_SFC_Scheafer = nanmean(brain_data_200p_Scheafer{:,2:201},1);

[r p] = corr(parcel_to_surface(averaged_SFC([181:360,1:180]),'glasser_360_fsa5')',parcel_to_surface(averaged_SFC_Scheafer,'schaefer_200_fsa5')')
X = zscore(averaged_SFC_Scheafer);[a, cb] = plot_cortical(parcel_to_surface(X,'schaefer_200_fsa5'),'color_range',[-2.7,2.7]);

figure;
scatter(parcel_to_surface(averaged_SFC([181:360,1:180]),'glasser_360_fsa5'), parcel_to_surface(averaged_SFC_Scheafer,'schaefer_200_fsa5'), 30, [0.5, 0.5, 0.5], 'filled', 'MarkerFaceAlpha', 0.3) % 半透明灰色
hold on;
p1 = polyfit(parcel_to_surface(averaged_SFC([181:360,1:180]),'glasser_360_fsa5')', parcel_to_surface(averaged_SFC_Scheafer,'schaefer_200_fsa5')', 1); % Linear fit
f1 = polyval(p1, parcel_to_surface(averaged_SFC([181:360,1:180]),'glasser_360_fsa5')');
plot(parcel_to_surface(averaged_SFC([181:360,1:180]),'glasser_360_fsa5')', f1, 'k-', 'LineWidth', 1); % Black fit line
ylim([0, 0.25]);
xlim([0, 0.25]);
hold off;


[corr_matrix, p_matrix] = corr(brain_data_7n{s2f,2:8}, brain_data_7n_Scheafer{f2s,2:8}, 'Rows', 'pairwise');
ICC = [diag(corr_matrix),diag(p_matrix)];
h3 = heatmap(Yeo7Names2, Yeo7Names2, corr_matrix, 'ColorbarVisible', 'on', 'Colormap', hot(100)); % Using heatmap function
h3.XLabel = 'SFC';
h3.YLabel = 'DSFC';
h3.Title = 'Correlation Heatmap';


[r p] = corr(parcel_to_surface(averaged_DSFC([181:360,1:180]),'glasser_360_fsa5')',parcel_to_surface(averaged_DSFC_Scheafer,'schaefer_200_fsa5')')
X = zscore(averaged_DSFC_Scheafer);[a, cb] = plot_cortical(parcel_to_surface(X,'schaefer_200_fsa5'),'color_range',[-2.7,2.7]);

figure;
scatter(parcel_to_surface(averaged_DSFC([181:360,1:180]),'glasser_360_fsa5'), parcel_to_surface(averaged_DSFC_Scheafer,'schaefer_200_fsa5'), 30, [0.5, 0.5, 0.5], 'filled', 'MarkerFaceAlpha', 0.3) % 半透明灰色
hold on;
p1 = polyfit(parcel_to_surface(averaged_DSFC([181:360,1:180]),'glasser_360_fsa5')', parcel_to_surface(averaged_DSFC_Scheafer,'schaefer_200_fsa5')', 1); % Linear fit
f1 = polyval(p1, parcel_to_surface(averaged_DSFC([181:360,1:180]),'glasser_360_fsa5')');
plot(parcel_to_surface(averaged_DSFC([181:360,1:180]),'glasser_360_fsa5')', f1, 'k-', 'LineWidth', 1); % Black fit line
ylim([0.07, 0.14]);
xlim([0.06, 0.12]);
hold off;

[corr_matrix, p_matrix] = corr(brain_data_7n{s2f,9:15}, brain_data_7n_Scheafer{f2s,9:15}, 'Rows', 'pairwise');
ICC = [diag(corr_matrix),diag(p_matrix)];
h3 = heatmap(Yeo7Names2, Yeo7Names2, corr_matrix, 'ColorbarVisible', 'on', 'Colormap', hot(100)); % Using heatmap function
h3.XLabel = 'SFC';
h3.YLabel = 'DSFC';
h3.Title = 'Correlation Heatmap';

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/ukb_glasser_scfcc_2_0_mean_FA.mat', 'SCDFCC_Mean', 'SCDFCC_SD', 'common_eid')
SFC_7n_FA = netMean1D(SCDFCC_Mean,Yeo7MMP);
DSFC_7n_FA = netMean1D(SCDFCC_SD,Yeo7MMP);

brain_data_7n_FA = [array2table(common_eid),array2table(SFC_7n_FA),array2table(DSFC_7n_FA)];
brain_data_7n_FA.Properties.VariableNames{1} = 'eid';
brain_data_360p_FA = [array2table(common_eid),array2table(SCDFCC_Mean),array2table(SCDFCC_SD)];
brain_data_360p_FA.Properties.VariableNames{1} = 'eid';

averaged_DSFC = nanmean(brain_data_360p{:,362:721},1);
averaged_SFC = nanmean(brain_data_360p{:,2:361},1);
scatter(averaged_SFC, averaged_DSFC, 5, Yeo7MMP, 'filled')
scatter(averaged_SFC, averaged_DSFC, 20, Yeo7MMP, 'filled')
Yeo7Colormap = [0.471,0.0710,0.522;0.275,0.510,0.706;0,0.463,0.0550;...
    0.769,0.224,0.976;0.863,0.973,0.639;0.902,0.576,0.129;0.804,0.239,0.306];
colormap(Yeo7Colormap)
scatter(averaged_SFC, averaged_DSFC, 50, Yeo7MMP, 'filled')
colormap(Yeo7Colormap)
scatter(netMean1D(averaged_SFC,Yeo7MMP), netMean1D(averaged_DSFC,Yeo7MMP), 500, 1:7, 'filled')
xlim([0, 0.25]);
ylim([0.06, 0.12]);
Corr_matrix = corr(SCDFCC_Mean, SCDFCC_SD,'Rows','pairwise');


Yeo7Names2 = {'VN','SMN','DAN','VAN','LN','FPN','DMN'};
Yeo7Colormap = [0.471,0.0710,0.522;0.275,0.510,0.706;0,0.463,0.0550;...
    0.769,0.224,0.976;0.863,0.973,0.639;0.902,0.576,0.129;0.804,0.239,0.306];

averaged_DSFC_FA = nanmean(brain_data_360p_FA{:,362:721},1);
averaged_SFC_FA = nanmean(brain_data_360p_FA{:,2:361},1);


parfor i = 1:size(brain_data_360p,1)
    [r_SFC(i,1),p_SFC(i,1)] = corr(brain_data_360p{i,2:361}',brain_data_360p_FA{i,2:361}');
    [r_DSFC(i,1),p_DSFC(i,1)] = corr(brain_data_360p{i,362:721}',brain_data_360p_FA{i,362:721}');
    disp(num2str(i))
end


figure;

% Subplot 1: Scatter plot with colormap and black fit line
subplot(2,3,1);
scatter(averaged_SFC, averaged_SFC_FA, 10, Yeo7MMP, 'filled');
colormap(Yeo7Colormap);
title('SFC');
hold on;
p1 = polyfit(averaged_SFC, averaged_SFC_FA, 1); % Linear fit
f1 = polyval(p1, averaged_SFC);
plot(averaged_SFC, f1, 'k-', 'LineWidth', 1); % Black fit line
hold off;

% Subplot 2: Histogram with black mean vertical line
subplot(2,3,2);
histogram(r_SFC, 50);
title('SFC');
hold on;
mean_val_SFC = nanmean(r_SFC);
xline(mean_val_SFC, 'k', 'LineWidth', 1); % Black vertical line at the mean
hold off;

% Subplot 4: Scatter plot with colormap and black fit line
subplot(2,3,4);
scatter(averaged_DSFC, averaged_DSFC_FA, 10, Yeo7MMP, 'filled');
colormap(Yeo7Colormap);
hold on;
p2 = polyfit(averaged_DSFC, averaged_DSFC_FA, 1); % Linear fit
f2 = polyval(p2, averaged_DSFC);
plot(averaged_DSFC, f2, 'k-', 'LineWidth', 1); % Black fit line
hold off;
title('DSFC');

% Subplot 5: Histogram with black mean vertical line
subplot(2,3,5);
histogram(r_DSFC, 50);
title('DSFC');
hold on;
mean_val_DSFC = nanmean(r_DSFC);
xline(mean_val_DSFC, 'k', 'LineWidth', 1); % Black vertical line at the mean
hold off;

% Subplot 3: Heatmap with custom colormap (using heatmap function)
subplot(2,3,3);
[corr_matrix, p_matrix] = corr(brain_data_7n{:,2:8}, brain_data_7n_FA{:,2:8}, 'Rows', 'pairwise');
ICC = [diag(corr_matrix),diag(p_matrix)];
h3 = heatmap(Yeo7Names2, Yeo7Names2, corr_matrix, 'ColorbarVisible', 'on', 'Colormap', hot(100)); % Using heatmap function
h3.XLabel = 'SFC';
h3.YLabel = 'DSFC';
h3.Title = 'Correlation Heatmap';

% Subplot 6: Heatmap with custom colormap (using heatmap function)
subplot(2,3,6);
[corr_matrix_2, p_matrix_2] = corr(brain_data_7n{:,9:15}, brain_data_7n_FA{:,9:15}, 'Rows', 'pairwise');
ICC = [diag(corr_matrix_2),diag(p_matrix_2)];
h6 = heatmap(Yeo7Names2, Yeo7Names2, corr_matrix_2, 'ColorbarVisible', 'on', 'Colormap', hot(100)); % Using heatmap function
h6.XLabel = 'SFC';
h6.YLabel = 'DSFC';
h6.Title = 'Correlation Heatmap';

[r p] = corr(averaged_SFC', averaged_SFC_FA')
[r p] = corr(averaged_DSFC', averaged_DSFC_FA')

%% Averaged Pattern Decoding
CytoColor = [110,45,125;73,104,153;180,210,143;240,204,49;240,233,61]./255;
Yeo7Colormap = [0.471,0.0710,0.522;0.275,0.510,0.706;0,0.463,0.0550;...
    0.769,0.224,0.976;0.863,0.973,0.639;0.902,0.576,0.129;0.804,0.239,0.306];

CytoName = {'Agranular','Frontal','Parietal','Polar','Granular'};
Yeo7Names2 = {'VN','SMN','DAN','VAN','LN','FPN','DMN'};

% Cytoarchitectonics 
%  i) agranular (purple; thick with large cells but sparse layers II and IV) 
%  ii) frontal (blue; thick but not rich in cellular substance, visible layers II and IV)
%  iii) parietal (green; thick and rich in cells with dense layers II and IV but small and slender pyramidal cells)
%  iv) polar (orange; thin but rich in cells, particularly in granular layers)
%  v) granular or koniocortex (yellow; thin but rich in smalls cells, even in layer IV, and a rarified layer V).
CytoAtlas   = dlmread(['economo_koskinas_fsa5.csv']);
GlasserAtlas   = dlmread(['glasser_360_fsa5.csv']);

for i = 1:360
    for j = 1:6
        dice_matrix(i,j) = dice(GlasserAtlas==i,CytoAtlas==j);
    end
end
for i = 1:360
    MMP2Cyto(i,1) = find(dice_matrix(i,:) == max(dice_matrix(i,:)));
end

figure;
subplot(2,2,1);b1 = bar(netMean1D(averaged_SFC,Yeo7MMP));
xticks(1:7); xticklabels(Yeo7Names2); xtickangle(90);ylim([.03,.15]);
b1.FaceColor = 'flat'; 
b1.CData = Yeo7Colormap;
subplot(2,2,2);b1 = bar(netMean1D(averaged_SFC,MMP2Cyto));
xticks(1:5); xticklabels(CytoName); xtickangle(90);ylim([.03,.15]);
b1.FaceColor = 'flat'; 
b1.CData = CytoColor;
subplot(2,2,3);b1 = bar(netMean1D(averaged_DSFC,Yeo7MMP));
xticks(1:7); xticklabels(Yeo7Names2); xtickangle(90);ylim([.05,.1]);
b1.FaceColor = 'flat'; 
b1.CData = Yeo7Colormap;
subplot(2,2,4);b1 = bar(netMean1D(averaged_DSFC,MMP2Cyto));
xticks(1:5); xticklabels(CytoName); xtickangle(90);ylim([.05,.1]);
b1.FaceColor = 'flat'; 
b1.CData = CytoColor;

Yeo17Names = {'VisualA';'VisualB';'SomatomotorA';'SomatomotorB';'DorsalAttentionA';'DorsalAttentionB';'SalienceA';...
    'SalienceB';'LimbicA';'LimbicB';'ControlC';'ControlA';'ControlB';'TemporalParietal';'DefaultC';'DefaultA';...
    'DefaultB'};

figure;
subplot(2,1,1);b1 = bar(netMean1D(averaged_SFC,Yeo17MMP));
xticks(1:17); xticklabels(Yeo17Names); xtickangle(90);ylim([0,.17]);
b1.FaceColor = 'flat'; 
b1.CData = rsn17cmap(2:18,:);
subplot(2,1,2);b1 = bar(netMean1D(averaged_DSFC,Yeo17MMP));
xticks(1:17); xticklabels(Yeo17Names); xtickangle(90);ylim([.05,.11]);
b1.FaceColor = 'flat'; 
b1.CData = rsn17cmap(2:18,:);

Yeo17Names = {'VisualA';'VisualB';'SomatomotorA';'SomatomotorB';'DorsalAttentionA';'DorsalAttentionB';'SalienceA';...
    'SalienceB';'LimbicA';'LimbicB';'ControlC';'ControlA';'ControlB';'TemporalParietal';'DefaultC';'DefaultA';...
    'DefaultB'};

figure;
subplot(2,1,1);b1 = bar(netMean1D(averaged_SFC,cole12_netkey));
xticks(1:12); xticklabels(cole12_name); xtickangle(90);ylim([0,.15]);
b1.FaceColor = 'flat'; 
b1.CData = cole12_color;
subplot(2,1,2);b1 = bar(netMean1D(averaged_DSFC,cole12_netkey));
xticks(1:12); xticklabels(cole12_name); xtickangle(90);ylim([.07,.1]);
b1.FaceColor = 'flat'; 
b1.CData = cole12_color;


% 360 regions
UKB_BrainXY = innerjoin(innerjoin(UKB_Basic,UKB_BrainMRICov),brain_data_360p);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_360p.Properties.VariableNames(2:end)'; 

formula = ['AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+t1_vol_scaling_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0'};
results_age_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
UKB_BrainXY.AgeAttend_2_0_Sq = zscore(UKB_BrainXY.AgeAttend_2_0).^2;
formula = ['AgeAttend_2_0+Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+t1_vol_scaling_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0_Sq'};
results_age2_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+SNR_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+t1_vol_scaling_2_0+(1|Centre_0_0)'];
XList = {'Sex_0_0'};
results_sex_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+HeadMotion_2_0+SNR_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+t1_vol_scaling_2_0+(1|Centre_0_0)'];
XList = {'BMI_2_0'};
results_bmi_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+SNR_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+(1|Centre_0_0)'];
XList = {'t1_vol_scaling_2_0'};
results_tiv_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
results_demo_360p = [results_age_360p;results_age2_360p;results_sex_360p;results_bmi_360p;results_tiv_360p];







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


brain_data_7n_Scheafer = [array2table(common_eid),array2table(SFC_7n_Scheafer),array2table(DSFC_7n_Scheafer)];
brain_data_7n_Scheafer.Properties.VariableNames{1} = 'eid';
brain_data_200p_Scheafer = [array2table(common_eid),array2table(SCDFCC_Mean),array2table(SCDFCC_SD)];
brain_data_200p_Scheafer.Properties.VariableNames{1} = 'eid';

[s2f,f2s] = idFinderNum(brain_data_360p.eid,brain_data_200p_Scheafer.eid);

[r_7n, p_7n] = corr(brain_data_7n{s2f,2:15}, brain_data_7n_Scheafer{f2s,2:15},'Rows','pairwise');


load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Structure_Function_Coupling_Aging/ukb_glasser_scfcc_2_0_mean_FA.mat', 'SCDFCC_Mean', 'SCDFCC_SD', 'common_eid')
SFC_7n_FA = netMean1D(SCDFCC_Mean,Yeo7MMP);
DSFC_7n_FA = netMean1D(SCDFCC_SD,Yeo7MMP);

brain_data_7n_FA = [array2table(common_eid),array2table(SFC_7n_FA),array2table(DSFC_7n_FA)];
brain_data_7n_FA.Properties.VariableNames{1} = 'eid';
brain_data_360p_FA = [array2table(common_eid),array2table(SCDFCC_Mean),array2table(SCDFCC_SD)];
brain_data_360p_FA.Properties.VariableNames{1} = 'eid';


C = [2:361]+360;
X1 = nanmean(brain_data_360p{:,C})';
X2 = nanmean(brain_data_360p_40TR{:,C})';
r = corr(X1,X2,'Rows','pairwise','Type','spearman')
p_perm = perm_sphere_p(X1,X2,perm_id,'spearman')

X1 = nanmean(brain_data_360p{:,C})';
X2 = nanmean(brain_data_360p_50TR{:,C})';
r = corr(X1,X2,'Rows','pairwise','Type','spearman')
p_perm = perm_sphere_p(X1,X2,perm_id,'spearman')

X1 = nanmean(brain_data_360p{:,C})';
X2 = nanmean(brain_data_360p_FA{:,C})';
r = corr(X1,X2,'Rows','pairwise','Type','spearman')
p_perm = perm_sphere_p(X1,X2,perm_id,'spearman')

X1 = nanmean(brain_data_360p{:,[2:361]+360})';
X2 = nanmean(brain_data_360p_50TR{:,2:361})';
r = corr(X1,X2,'Rows','pairwise','Type','spearman')
p_perm = perm_sphere_p(X1,X2,perm_id,'spearman')


