function results = run_genetic_associations(window_length)
% RUN_GENETIC_ASSOCIATIONS
%
% Reproduce PRS, APOE dosage, AD-genotype and age-by-genetic-risk analyses
% from the original MATLAB master script.
%
% Required tables:
%   UKB_Basic.csv
%   UKB_BrainMRICov.csv
%   UKB_GeneCov.csv
%   UKB_PRS.tsv
%   APOE4_genetype.txt
%   AD_genetype.txt
%   brain_data_<window>TR.mat
%
% No participant-level genetic data are distributed with this repository.
%
% -------------------------------------------------------------------------

if nargin < 1 || isempty(window_length)
    window_length = 30;
end

run('00_config.m');

table_dir = fullfile(OUTPUT_DIR, 'analysis_tables');

brain_file = fullfile( ...
    table_dir, ...
    sprintf('brain_data_%dTR.mat', window_length));

load(brain_file, 'brain_data_7n', 'brain_data_360p');

UKB_Basic = readtable(fullfile(DATA_ROOT, 'tabular', 'UKB_Basic.csv'));
UKB_BrainMRICov = readtable(fullfile(DATA_ROOT, 'tabular', 'UKB_BrainMRICov.csv'));

UKB_GeneCov = readtable( ...
    fullfile(DATA_ROOT, 'genetics', 'UKB_GeneCov.csv'));

UKB_PRS = readtable( ...
    fullfile(DATA_ROOT, 'genetics', 'UKB_PRS.tsv'), ...
    'FileType', ...
    'text');

APOE4_genetype = readtable( ...
    fullfile(DATA_ROOT, 'genetics', 'APOE4_genetype.txt'));

AD_genetype = readtable( ...
    fullfile(DATA_ROOT, 'genetics', 'AD_genetype.txt'));

% Harmonize participant-ID variable names where possible.
UKB_GeneCov = harmonize_eid(UKB_GeneCov);
UKB_PRS = harmonize_eid(UKB_PRS);
APOE4_genetype = harmonize_eid(APOE4_genetype);
AD_genetype = harmonize_eid(AD_genetype);

results = struct();

genetic_formula = [ ...
    'AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+', ...
    'HeadMotion_2_0+SNR_2_0+t1_vol_scaling_2_0+', ...
    'x22000_0_0+', ...
    'x22009_0_1+x22009_0_2+x22009_0_3+x22009_0_4+', ...
    'x22009_0_5+x22009_0_6+x22009_0_7+x22009_0_8+', ...
    'x22009_0_9+x22009_0_10+', ...
    '(1|Centre_0_0)' ...
];

% -------------------------------------------------------------------------
% Yeo-7 PRS
% -------------------------------------------------------------------------

UKB_BrainXY = innerjoin(UKB_Basic, UKB_BrainMRICov, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, UKB_GeneCov, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, UKB_PRS, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, brain_data_7n, 'Keys', 'eid');

UKB_BrainXY = apply_motion_filter(UKB_BrainXY);

YList = brain_data_7n.Properties.VariableNames(2:end)';
XList = setdiff(UKB_PRS.Properties.VariableNames, {'eid'}, 'stable')';

results.PRS_7n = LMM_with_EffectSizes_Full( ...
    UKB_BrainXY, ...
    genetic_formula, ...
    XList, ...
    YList);

% -------------------------------------------------------------------------
% Yeo-7 APOE dosage
% -------------------------------------------------------------------------

UKB_BrainXY = innerjoin(UKB_Basic, UKB_BrainMRICov, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, UKB_GeneCov, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, APOE4_genetype, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, brain_data_7n, 'Keys', 'eid');

UKB_BrainXY = apply_motion_filter(UKB_BrainXY);

XList = {'APOE4_dosage'};

results.APOE4_7n = LMM_with_EffectSizes_Full( ...
    UKB_BrainXY, ...
    genetic_formula, ...
    XList, ...
    YList);

% -------------------------------------------------------------------------
% Yeo-7 AD genotype variables
% -------------------------------------------------------------------------

UKB_BrainXY_AD = innerjoin(UKB_Basic, UKB_BrainMRICov, 'Keys', 'eid');
UKB_BrainXY_AD = innerjoin(UKB_BrainXY_AD, UKB_GeneCov, 'Keys', 'eid');
UKB_BrainXY_AD = innerjoin(UKB_BrainXY_AD, AD_genetype, 'Keys', 'eid');
UKB_BrainXY_AD = innerjoin(UKB_BrainXY_AD, brain_data_7n, 'Keys', 'eid');

UKB_BrainXY_AD = apply_motion_filter(UKB_BrainXY_AD);

XList_AD = setdiff(AD_genetype.Properties.VariableNames, {'eid'}, 'stable')';

results.AD_genes_7n = LMM_with_EffectSizes_Full( ...
    UKB_BrainXY_AD, ...
    genetic_formula, ...
    XList_AD, ...
    YList);

% -------------------------------------------------------------------------
% APOE dosage x age interaction
% -------------------------------------------------------------------------

interaction_formula = [ ...
    'Sex_0_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+SNR_2_0+', ...
    't1_vol_scaling_2_0+x22000_0_0+', ...
    'x22009_0_1+x22009_0_2+x22009_0_3+x22009_0_4+', ...
    'x22009_0_5+x22009_0_6+x22009_0_7+x22009_0_8+', ...
    'x22009_0_9+x22009_0_10+', ...
    '(1|Centre_0_0)' ...
];

X1List = {'APOE4_dosage'};
X2List = {'AgeAttend_2_0'};

results.APOE_age_interaction_7n = LMM_withFormula_Interact( ...
    UKB_BrainXY, ...
    interaction_formula, ...
    X1List, ...
    X2List, ...
    YList);

% -------------------------------------------------------------------------
% Parcel-level PRS
% -------------------------------------------------------------------------

UKB_BrainXY = innerjoin(UKB_Basic, UKB_BrainMRICov, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, UKB_GeneCov, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, UKB_PRS, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, brain_data_360p, 'Keys', 'eid');

UKB_BrainXY = apply_motion_filter(UKB_BrainXY);

YList360 = brain_data_360p.Properties.VariableNames(2:end)';

% Preserve the two PRS variables explicitly used by the original parcel-level
% analysis when present.
preferred_parcel_prs = {'p26206', 'p26275'};
XList360 = preferred_parcel_prs( ...
    ismember(preferred_parcel_prs, ...
    UKB_BrainXY.Properties.VariableNames));

if isempty(XList360)
    warning(['p26206/p26275 were not found. Parcel-level PRS analysis ', ...
             'will be skipped rather than substituting other PRS variables.']);
    results.PRS_360p = [];
else
    results.PRS_360p = LMM_with_EffectSizes_Full( ...
        UKB_BrainXY, ...
        genetic_formula, ...
        XList360, ...
        YList360);
end

% -------------------------------------------------------------------------
% Parcel-level APOE
% -------------------------------------------------------------------------

UKB_BrainXY = innerjoin(UKB_Basic, UKB_BrainMRICov, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, UKB_GeneCov, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, APOE4_genetype, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, brain_data_360p, 'Keys', 'eid');

UKB_BrainXY = apply_motion_filter(UKB_BrainXY);

results.APOE4_360p = LMM_with_EffectSizes_Full( ...
    UKB_BrainXY, ...
    genetic_formula, ...
    {'APOE4_dosage'}, ...
    YList360);

results.APOE_age_interaction_360p = LMM_withFormula_Interact( ...
    UKB_BrainXY, ...
    interaction_formula, ...
    {'APOE4_dosage'}, ...
    {'AgeAttend_2_0'}, ...
    YList360);

% -------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------

result_dir = fullfile(OUTPUT_DIR, 'association_results');

if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

save( ...
    fullfile(result_dir, ...
    sprintf('genetic_results_%dTR.mat', window_length)), ...
    'results', ...
    '-v7.3');

end


function T = harmonize_eid(T)
% Convert common historical participant-ID variable names to 'eid'.

vars = T.Properties.VariableNames;

candidates = {'eid', 'IID', 'iid', 'ID', 'participant_id'};

for i = 1:numel(candidates)
    if ismember(candidates{i}, vars)
        if ~strcmp(candidates{i}, 'eid')
            T.Properties.VariableNames{strcmp(vars, candidates{i})} = 'eid';
        end
        return;
    end
end

error('No participant-ID column could be identified.');
end


function T = apply_motion_filter(T)

if ismember('HeadMotion_2_0', T.Properties.VariableNames)
    T(T.HeadMotion_2_0 > 0.2, :) = [];
end

end
