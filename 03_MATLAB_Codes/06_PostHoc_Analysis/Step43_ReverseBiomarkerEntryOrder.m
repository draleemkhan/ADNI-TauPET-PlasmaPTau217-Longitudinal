%% ============================================================
% STEP 43: REVERSE BIOMARKER ENTRY-ORDER ANALYSIS
%
% Purpose:
%   Determine whether the incremental contribution of tau PET
%   and plasma p-tau217 depends on biomarker entry order.
%
% Original sequence (Step 20):
%   Model 1: Centiloid
%   Model 2: + plasma p-tau217
%   Model 3: + tau PET
%
% Reverse sequence (Step 43):
%   Model 1: Centiloid
%   Model 2R: + tau PET
%   Model 3: + plasma p-tau217
%
% All models:
%   - Maximum likelihood (ML)
%   - Same 379 observations / 143 participants
%   - Random intercept + uncorrelated random time slope
%   - Same demographic adjustment
%
% Outputs:
%   Step43_ReverseEntry_ModelComparison.csv
%   Step43_ReverseEntry_IncrementalTests.csv
%   Step43_ReverseEntry_Models.mat
% ============================================================

clear;
clc;

fprintf('\n============================================================\n');
fprintf(' STEP 43: REVERSE BIOMARKER ENTRY ORDER\n');
fprintf('============================================================\n');

%% ------------------------------------------------------------
% Paths
% -------------------------------------------------------------

outputDir = ...
    'D:\output';

dataFile = fullfile(outputDir,...
    'Step08_ModelReady_379_Observations.csv');

step20File = fullfile(outputDir,...
    'Step20_Incremental_Biomarker_Models_ML_Summary.csv');

%% ------------------------------------------------------------
% Load model-ready data
% -------------------------------------------------------------

T = readtable(dataFile,...
    'VariableNamingRule','preserve');

fprintf('\nDataset: %d observations\n',height(T));

%% ------------------------------------------------------------
% Required-variable audit
% -------------------------------------------------------------

required = [...
    "RID"
    "CDRSB"
    "TimeYears"
    "Age"
    "Education"
    "Tau_z"
    "Centiloid_z"
    "PTau217_z"
    "Sex"];

for k = 1:numel(required)
    assert(ismember(required(k),string(T.Properties.VariableNames)),...
        'Missing required variable: %s',required(k));
end

nParticipants = numel(unique(T.RID));

fprintf('Participants: %d\n',nParticipants);

assert(height(T)==379,...
    'Expected 379 observations.');

assert(nParticipants==143,...
    'Expected 143 participants.');

fprintf('Dataset audit: PASS\n');

%% ------------------------------------------------------------
% Reconstruct sex coding exactly as primary model
% Male = reference
% -------------------------------------------------------------

T.SexModel = categorical(string(T.Sex));

cats = categories(T.SexModel);

assert(any(strcmp(cats,'Male')),...
    'Male category not found.');

assert(any(strcmp(cats,'Female')),...
    'Female category not found.');

T.SexModel = reordercats(T.SexModel,...
    {'Male','Female'});

fprintf('\nSex reference category: Male\n');

%% ------------------------------------------------------------
% Model formulas
%
% Important:
% Each biomarker enters with BOTH its main effect and
% biomarker x time interaction.
%
% Random structure:
% (1 + TimeYears | RID) with covariance pattern diagonal
% implemented as:
% (1 | RID) + (TimeYears - 1 | RID)
% -------------------------------------------------------------

formulaM1 = [...
    'CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
    'Centiloid_z + TimeYears:Centiloid_z + ' ...
    '(1|RID) + (TimeYears-1|RID)'];

formulaM2R = [...
    'CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
    'Centiloid_z + TimeYears:Centiloid_z + ' ...
    'Tau_z + TimeYears:Tau_z + ' ...
    '(1|RID) + (TimeYears-1|RID)'];

formulaM3 = [...
    'CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
    'Centiloid_z + TimeYears:Centiloid_z + ' ...
    'Tau_z + TimeYears:Tau_z + ' ...
    'PTau217_z + TimeYears:PTau217_z + ' ...
    '(1|RID) + (TimeYears-1|RID)'];

%% ------------------------------------------------------------
% Fit models using maximum likelihood
% -------------------------------------------------------------

fprintf('\nFitting Model 1: Centiloid...\n');

mdl1 = fitlme(T,formulaM1,...
    'FitMethod','ML');

fprintf('Fitting Model 2R: + tau PET...\n');

mdl2R = fitlme(T,formulaM2R,...
    'FitMethod','ML');

fprintf('Fitting Model 3: + plasma p-tau217...\n');

mdl3 = fitlme(T,formulaM3,...
    'FitMethod','ML');

%% ------------------------------------------------------------
% Extract model criteria
% -------------------------------------------------------------

IC1  = mdl1.ModelCriterion;
IC2R = mdl2R.ModelCriterion;
IC3  = mdl3.ModelCriterion;

LL1  = mdl1.LogLikelihood;
LL2R = mdl2R.LogLikelihood;
LL3  = mdl3.LogLikelihood;

%% ------------------------------------------------------------
% Likelihood-ratio tests
%
% M1 -> M2R:
% Incremental contribution of tau PET after Centiloid
%
% M2R -> M3:
% Incremental contribution of p-tau217 after tau PET
% -------------------------------------------------------------

cmpTau = compare(mdl1,mdl2R);
cmpPTau = compare(mdl2R,mdl3);

disp(' ');
disp('============================================================');
disp(' ADD TAU PET AFTER CENTILOID');
disp('============================================================');
disp(cmpTau);

disp(' ');
disp('============================================================');
disp(' ADD PLASMA p-tau217 AFTER TAU PET');
disp('============================================================');
disp(cmpPTau);

%% ------------------------------------------------------------
% Calculate LR tests directly for transparent output
% Each biomarker block adds two parameters:
% main effect + biomarker x time interaction
% -------------------------------------------------------------

LR_Tau = 2*(LL2R-LL1);
df_Tau = mdl2R.NumCoefficients - mdl1.NumCoefficients;
p_Tau  = 1-chi2cdf(LR_Tau,df_Tau);

LR_PTau = 2*(LL3-LL2R);
df_PTau = mdl3.NumCoefficients - mdl2R.NumCoefficients;
p_PTau  = 1-chi2cdf(LR_PTau,df_PTau);

%% ------------------------------------------------------------
% Changes in information criteria
% Negative values favour the larger model
% -------------------------------------------------------------

DeltaAIC_Tau = IC2R.AIC - IC1.AIC;
DeltaBIC_Tau = IC2R.BIC - IC1.BIC;

DeltaAIC_PTau = IC3.AIC - IC2R.AIC;
DeltaBIC_PTau = IC3.BIC - IC2R.BIC;

%% ------------------------------------------------------------
% Model-comparison table
% -------------------------------------------------------------

Model = [...
    "Model 1"
    "Model 2R"
    "Model 3"];

BiomarkerBlock = [...
    "Centiloid"
    "Centiloid + tau PET"
    "Centiloid + tau PET + plasma p-tau217"];

AIC = [...
    IC1.AIC
    IC2R.AIC
    IC3.AIC];

BIC = [...
    IC1.BIC
    IC2R.BIC
    IC3.BIC];

LogLikelihood = [...
    LL1
    LL2R
    LL3];

ModelComparison = table(...
    Model,...
    BiomarkerBlock,...
    AIC,...
    BIC,...
    LogLikelihood);

%% ------------------------------------------------------------
% Incremental-test table
% -------------------------------------------------------------

Addition = [...
    "Tau PET after Centiloid"
    "Plasma p-tau217 after Centiloid + tau PET"];

LRStat = [...
    LR_Tau
    LR_PTau];

DF = [...
    df_Tau
    df_PTau];

PValue = [...
    p_Tau
    p_PTau];

DeltaAIC = [...
    DeltaAIC_Tau
    DeltaAIC_PTau];

DeltaBIC = [...
    DeltaBIC_Tau
    DeltaBIC_PTau];

IncrementalTests = table(...
    Addition,...
    LRStat,...
    DF,...
    PValue,...
    DeltaAIC,...
    DeltaBIC);

%% ------------------------------------------------------------
% Display results
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' REVERSE-SEQUENCE MODEL FIT\n');
fprintf('============================================================\n\n');

disp(ModelComparison);

fprintf('\n============================================================\n');
fprintf(' INCREMENTAL TESTS\n');
fprintf('============================================================\n\n');

disp(IncrementalTests);

fprintf('\nTau PET added after Centiloid:\n');
fprintf('  LR chi-square(%d) = %.4f\n',df_Tau,LR_Tau);
fprintf('  P = %.8g\n',p_Tau);
fprintf('  Delta AIC = %.4f\n',DeltaAIC_Tau);
fprintf('  Delta BIC = %.4f\n',DeltaBIC_Tau);

fprintf('\nPlasma p-tau217 added after Tau PET:\n');
fprintf('  LR chi-square(%d) = %.4f\n',df_PTau,LR_PTau);
fprintf('  P = %.8g\n',p_PTau);
fprintf('  Delta AIC = %.4f\n',DeltaAIC_PTau);
fprintf('  Delta BIC = %.4f\n',DeltaBIC_PTau);

%% ------------------------------------------------------------
% Cross-check final full model against authoritative Step 20
% -------------------------------------------------------------

S20 = readtable(step20File,...
    'VariableNamingRule','preserve');

assert(height(S20)==1,...
    'Unexpected Step 20 summary structure.');

authAIC3 = S20.Model3_AIC(1);
authBIC3 = S20.Model3_BIC(1);

diffAIC = abs(IC3.AIC-authAIC3);
diffBIC = abs(IC3.BIC-authBIC3);

fprintf('\n============================================================\n');
fprintf(' FULL-MODEL REPRODUCTION CHECK\n');
fprintf('============================================================\n');

fprintf('Step 20 Model 3 AIC: %.6f\n',authAIC3);
fprintf('Step 43 Model 3 AIC: %.6f\n',IC3.AIC);
fprintf('Absolute difference: %.10f\n',diffAIC);

fprintf('\nStep 20 Model 3 BIC: %.6f\n',authBIC3);
fprintf('Step 43 Model 3 BIC: %.6f\n',IC3.BIC);
fprintf('Absolute difference: %.10f\n',diffBIC);

% Step20 CSV values are rounded, so allow small tolerance
assert(diffAIC < 0.1,...
    'Final-model AIC does not reproduce Step 20.');

assert(diffBIC < 0.1,...
    'Final-model BIC does not reproduce Step 20.');

fprintf('\nFINAL FULL MODEL REPRODUCTION: PASS\n');

%% ------------------------------------------------------------
% Save outputs
% -------------------------------------------------------------

writetable(ModelComparison,...
    fullfile(outputDir,...
    'Step43_ReverseEntry_ModelComparison.csv'));

writetable(IncrementalTests,...
    fullfile(outputDir,...
    'Step43_ReverseEntry_IncrementalTests.csv'));

save(fullfile(outputDir,...
    'Step43_ReverseEntry_Models.mat'),...
    'mdl1','mdl2R','mdl3',...
    'ModelComparison','IncrementalTests',...
    'cmpTau','cmpPTau');

%% ------------------------------------------------------------
% Complete
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' STEP 43 COMPLETE\n');
fprintf('============================================================\n');

fprintf('Observations: 379: PASS\n');
fprintf('Participants: 143: PASS\n');
fprintf('Reverse sequence: Centiloid -> Tau PET -> p-tau217\n');
fprintf('Final full-model reproduction: PASS\n');

fprintf('\nSaved outputs:\n');
fprintf('  Step43_ReverseEntry_ModelComparison.csv\n');
fprintf('  Step43_ReverseEntry_IncrementalTests.csv\n');
fprintf('  Step43_ReverseEntry_Models.mat\n');

fprintf('\n============================================================\n');