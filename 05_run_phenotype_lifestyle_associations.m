function results = run_phenotype_lifestyle_associations(window_length)
% RUN_PHENOTYPE_LIFESTYLE_ASSOCIATIONS
%
% Reproduce phenotype, cognition/health, and lifestyle association analyses
% from the original MATLAB master script.
%
% Required tables:
%   UKB_Basic.csv
%   UKB_BrainMRICov.csv
%   UKB_Outcomes_2_0.csv
%   UKB_Lifestyle_0_0.csv
%   brain_data_<window>TR.mat
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
UKB_Outcomes_2_0 = readtable(fullfile(DATA_ROOT, 'tabular', 'UKB_Outcomes_2_0.csv'));
UKB_Lifestyle_0_0 = readtable(fullfile(DATA_ROOT, 'tabular', 'UKB_Lifestyle_0_0.csv'));

results = struct();

% -------------------------------------------------------------------------
% Outcome associations: Yeo-7
% -------------------------------------------------------------------------

UKB_BrainXY = innerjoin(UKB_Basic, UKB_BrainMRICov, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, UKB_Outcomes_2_0, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, brain_data_7n, 'Keys', 'eid');

if ismember('HeadMotion_2_0', UKB_BrainXY.Properties.VariableNames)
    UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0 > 0.2, :) = [];
end

YList = brain_data_7n.Properties.VariableNames(2:end)';
XList = UKB_Outcomes_2_0.Properties.VariableNames(2:end)';

formula = [ ...
    'AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+', ...
    'HeadMotion_2_0+SNR_2_0+', ...
    'Race_1+Race_2+Race_3+Race_4+', ...
    't1_vol_scaling_2_0+(1|Centre_0_0)' ...
];

results.pheno_7n = LMM_with_EffectSizes_Full( ...
    UKB_BrainXY, ...
    formula, ...
    XList, ...
    YList);

% -------------------------------------------------------------------------
% Lifestyle associations: Yeo-7
% -------------------------------------------------------------------------

UKB_BrainXY = innerjoin(UKB_Basic, UKB_BrainMRICov, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, UKB_Lifestyle_0_0, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, brain_data_7n, 'Keys', 'eid');

if ismember('HeadMotion_2_0', UKB_BrainXY.Properties.VariableNames)
    UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0 > 0.2, :) = [];
end

YList = brain_data_7n.Properties.VariableNames(2:end)';

% Prefer the exact downstream variables used by the manuscript code if they
% exist; otherwise use all lifestyle variables except eid.
preferred_lifestyle_vars = { ...
    'x738', ...
    'x20022', ...
    'x20116', ...
    'x20495', ...
    'x20497', ...
    'x22189', ...
    'x1160_L' ...
};

missing_lifestyle = setdiff( ...
    preferred_lifestyle_vars, ...
    UKB_BrainXY.Properties.VariableNames);

if ~isempty(missing_lifestyle)
    error( ...
        ['The analysis-ready lifestyle table is missing variables used by ', ...
         'the original manuscript analysis: ', strjoin(missing_lifestyle, ', ')]);
end

XList = preferred_lifestyle_vars;

formula = [ ...
    'AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+', ...
    'HeadMotion_2_0+SNR_2_0+', ...
    'Race_1+Race_2+Race_3+Race_4+', ...
    't1_vol_scaling_2_0+(1|Centre_0_0)' ...
];

results.lifestyle_7n = LMM_with_EffectSizes_Full( ...
    UKB_BrainXY, ...
    formula, ...
    XList, ...
    YList);

% -------------------------------------------------------------------------
% Outcome associations: parcel-level
% -------------------------------------------------------------------------

UKB_BrainXY = innerjoin(UKB_Basic, UKB_BrainMRICov, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, UKB_Outcomes_2_0, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, brain_data_360p, 'Keys', 'eid');

if ismember('HeadMotion_2_0', UKB_BrainXY.Properties.VariableNames)
    UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0 > 0.2, :) = [];
end

YList = brain_data_360p.Properties.VariableNames(2:end)';

preferred_outcome_vars = { ...
    'x135', 'x137', 'x1970', 'x2010', ...
    'x2178', 'x2188', 'x4282', 'x4526', ...
    'x4548', 'x6373', 'x20016', 'x20197' ...
};

missing_outcomes = setdiff( ...
    preferred_outcome_vars, ...
    UKB_BrainXY.Properties.VariableNames);

if ~isempty(missing_outcomes)
    error( ...
        ['UKB_Outcomes_2_0 is missing variables used by the original ', ...
         'parcel-level analysis: ', strjoin(missing_outcomes, ', ')]);
end

XList = preferred_outcome_vars;

results.pheno_360p = LMM_with_EffectSizes_Full( ...
    UKB_BrainXY, ...
    formula, ...
    XList, ...
    YList);

% -------------------------------------------------------------------------
% Lifestyle associations: parcel-level
% -------------------------------------------------------------------------

UKB_BrainXY = innerjoin(UKB_Basic, UKB_BrainMRICov, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, UKB_Lifestyle_0_0, 'Keys', 'eid');
UKB_BrainXY = innerjoin(UKB_BrainXY, brain_data_360p, 'Keys', 'eid');

if ismember('HeadMotion_2_0', UKB_BrainXY.Properties.VariableNames)
    UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0 > 0.2, :) = [];
end

YList = brain_data_360p.Properties.VariableNames(2:end)';

missing_lifestyle = setdiff( ...
    preferred_lifestyle_vars, ...
    UKB_BrainXY.Properties.VariableNames);

if ~isempty(missing_lifestyle)
    error( ...
        ['The analysis-ready lifestyle table is missing variables used by ', ...
         'the original manuscript analysis: ', strjoin(missing_lifestyle, ', ')]);
end

XList = preferred_lifestyle_vars;

results.lifestyle_360p = LMM_with_EffectSizes_Full( ...
    UKB_BrainXY, ...
    formula, ...
    XList, ...
    YList);

% -------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------

result_dir = fullfile(OUTPUT_DIR, 'association_results');

if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

save( ...
    fullfile(result_dir, ...
    sprintf('phenotype_lifestyle_results_%dTR.mat', window_length)), ...
    'results', ...
    '-v7.3');

end
