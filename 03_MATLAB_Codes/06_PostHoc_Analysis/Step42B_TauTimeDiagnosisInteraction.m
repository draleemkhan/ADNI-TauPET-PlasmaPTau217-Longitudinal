%% Step42B_TauTimeDiagnosisInteraction.m
% ============================================================
% FORMAL DIAGNOSIS HETEROGENEITY ANALYSIS
%
% Tests whether the baseline Tau PET x Time association differs
% according to baseline diagnostic stage.
%
% Reference diagnosis: MCI
%
% Full model includes:
%   Primary-model terms
%   + Diagnosis
%   + Diagnosis x Time
%   + Diagnosis x Tau
%   + Diagnosis x Tau x Time
%
% Comparison:
%   Reduced model = same model without Tau x Time x Diagnosis
%   Full model    = includes Tau x Time x Diagnosis
%
% ML is used for the likelihood-ratio comparison.
% ============================================================

clear; clc;

fprintf('\n============================================================\n');
fprintf(' STEP 42B: TAU x TIME x DIAGNOSIS INTERACTION\n');
fprintf('============================================================\n');

outputDir = ...
    'D:\output';

%% ------------------------------------------------------------
% Load model-ready data
% -------------------------------------------------------------

T = readtable(fullfile(outputDir,...
    'Step08_ModelReady_379_Observations.csv'),...
    'VariableNamingRule','preserve');

fprintf('\nDataset: %d observations / %d participants\n',...
    height(T),numel(unique(T.RID)));

%% ------------------------------------------------------------
% Sex coding exactly as primary model
% -------------------------------------------------------------

T.SexModel = categorical(string(T.Sex),...
    ["Male","Female"]);

T.SexModel = reordercats(T.SexModel,...
    {'Male','Female'});

%% ------------------------------------------------------------
% Diagnosis coding
% MCI deliberately selected as reference category
% -------------------------------------------------------------

T.DxModel = categorical(string(T.Baseline_Diagnosis),...
    ["MCI","Cognitively normal","Dementia"]);

T.DxModel = reordercats(T.DxModel,...
    {'MCI','Cognitively normal','Dementia'});

fprintf('\nDiagnosis reference category: MCI\n');

% Participant-level audit
[~,ia] = unique(T.RID,'stable');

dxParticipant = T.DxModel(ia);

fprintf('\nParticipant-level diagnosis counts:\n');
disp(countcats(dxParticipant));
disp(categories(dxParticipant));

assert(sum(dxParticipant=='MCI')==75);
assert(sum(dxParticipant=='Cognitively normal')==33);
assert(sum(dxParticipant=='Dementia')==35);

fprintf('Diagnosis-count audit: PASS\n');

%% ------------------------------------------------------------
% REDUCED MODEL
%
% Contains all lower-order terms but excludes only the
% diagnosis-specific modification of Tau x Time.
%
% ML is required for nested-model likelihood-ratio testing.
% -------------------------------------------------------------

formulaReduced = [ ...
    'CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
    'Tau_z + Centiloid_z + PTau217_z + DxModel + ' ...
    'TimeYears:Tau_z + TimeYears:Centiloid_z + ' ...
    'TimeYears:PTau217_z + ' ...
    'TimeYears:DxModel + Tau_z:DxModel + ' ...
    '(1|RID) + (TimeYears-1|RID)' ...
    ];

fprintf('\nFitting reduced ML model...\n');

mdlReduced = fitlme(T,formulaReduced,...
    'FitMethod','ML');

%% ------------------------------------------------------------
% FULL MODEL
%
% Adds Tau x Time x Diagnosis
% -------------------------------------------------------------

formulaFull = [ ...
    'CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
    'Tau_z + Centiloid_z + PTau217_z + DxModel + ' ...
    'TimeYears:Tau_z + TimeYears:Centiloid_z + ' ...
    'TimeYears:PTau217_z + ' ...
    'TimeYears:DxModel + Tau_z:DxModel + ' ...
    'TimeYears:Tau_z:DxModel + ' ...
    '(1|RID) + (TimeYears-1|RID)' ...
    ];

fprintf('Fitting full ML model...\n');

mdlFull = fitlme(T,formulaFull,...
    'FitMethod','ML');

%% ------------------------------------------------------------
% Fixed effects
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' FULL MODEL FIXED EFFECTS\n');
fprintf('============================================================\n');

C = mdlFull.Coefficients;
disp(C);

%% ------------------------------------------------------------
% Likelihood-ratio test of the TWO three-way interaction terms
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' GLOBAL TAU x TIME x DIAGNOSIS TEST\n');
fprintf('============================================================\n');

cmp = compare(mdlReduced,mdlFull);

disp(cmp);

%% ------------------------------------------------------------
% Extract likelihood ratio information robustly
% -------------------------------------------------------------

LLred  = mdlReduced.LogLikelihood;
LLfull = mdlFull.LogLikelihood;

LR = 2*(LLfull-LLred);

dfDifference = mdlFull.NumCoefficients - ...
               mdlReduced.NumCoefficients;

pGlobal = 1 - chi2cdf(LR,dfDifference);

fprintf('\nReduced log-likelihood: %.6f\n',LLred);
fprintf('Full log-likelihood:    %.6f\n',LLfull);
fprintf('LR statistic:           %.6f\n',LR);
fprintf('df difference:          %d\n',dfDifference);
fprintf('GLOBAL interaction P:   %.6g\n',pGlobal);

assert(dfDifference==2,...
    'Expected two diagnosis-specific three-way interaction terms.');

%% ------------------------------------------------------------
% Identify Tau x Time and three-way coefficient names
% -------------------------------------------------------------

names = string(C.Name);

fprintf('\nRelevant Tau / Time / Diagnosis coefficients:\n');

idxRelevant = ...
    contains(names,'Tau_z') & ...
    (contains(names,'TimeYears') | contains(names,'DxModel'));

disp(C(idxRelevant,:));

%% ------------------------------------------------------------
% Identify reference-group Tau x Time coefficient
% -------------------------------------------------------------

baseIdx = strcmp(names,'TimeYears:Tau_z') | ...
          strcmp(names,'Tau_z:TimeYears');

assert(sum(baseIdx)==1,...
    'Could not uniquely identify reference Tau x Time coefficient.');

b = C.Estimate;
V = mdlFull.CoefficientCovariance;

betaMCI = b(baseIdx);

%% ------------------------------------------------------------
% Find diagnosis-specific three-way terms
% -------------------------------------------------------------

threeIdx = ...
    contains(names,'TimeYears') & ...
    contains(names,'Tau_z') & ...
    contains(names,'DxModel');

threeNames = names(threeIdx);

fprintf('\nThree-way interaction coefficient names:\n');
disp(threeNames);

assert(sum(threeIdx)==2,...
    'Expected exactly two Tau x Time x Diagnosis coefficients.');

cn3 = threeIdx & contains(names,'Cognitively normal');
dem3 = threeIdx & contains(names,'Dementia');

assert(sum(cn3)==1,...
    'Could not identify CN three-way coefficient.');

assert(sum(dem3)==1,...
    'Could not identify dementia three-way coefficient.');

%% ------------------------------------------------------------
% Linear combinations:
%
% MCI:
%   Tau x Time
%
% CN:
%   Tau x Time + CN three-way interaction
%
% Dementia:
%   Tau x Time + Dementia three-way interaction
% -------------------------------------------------------------

k = numel(b);

L_MCI = zeros(1,k);
L_CN  = zeros(1,k);
L_DEM = zeros(1,k);

L_MCI(baseIdx) = 1;

L_CN(baseIdx) = 1;
L_CN(cn3)     = 1;

L_DEM(baseIdx) = 1;
L_DEM(dem3)    = 1;

%% ------------------------------------------------------------
% Helper calculations
% -------------------------------------------------------------

estimateCI = @(L) deal(...
    L*b,...
    sqrt(L*V*L'),...
    L*b - 1.96*sqrt(L*V*L'),...
    L*b + 1.96*sqrt(L*V*L'));

[bMCI,seMCI,loMCI,hiMCI] = estimateCI(L_MCI);
[bCN,seCN,loCN,hiCN]      = estimateCI(L_CN);
[bDEM,seDEM,loDEM,hiDEM]  = estimateCI(L_DEM);

% Wald P values
pMCI = 2*(1-normcdf(abs(bMCI/seMCI)));
pCN  = 2*(1-normcdf(abs(bCN/seCN)));
pDEM = 2*(1-normcdf(abs(bDEM/seDEM)));

%% ------------------------------------------------------------
% Stage-specific conditional Tau x Time estimates
% -------------------------------------------------------------

StageSpecific = table(...
    ["MCI";"Cognitively normal";"Dementia"],...
    [75;33;35],...
    [bMCI;bCN;bDEM],...
    [seMCI;seCN;seDEM],...
    [loMCI;loCN;loDEM],...
    [hiMCI;hiCN;hiDEM],...
    [pMCI;pCN;pDEM],...
    'VariableNames',{...
    'Diagnosis','N','Beta','SE',...
    'CI_Lower','CI_Upper','PValue'});

fprintf('\n============================================================\n');
fprintf(' DIAGNOSIS-SPECIFIC TAU PET x TIME ESTIMATES\n');
fprintf(' FROM FORMAL INTERACTION MODEL\n');
fprintf('============================================================\n');

disp(StageSpecific);

%% ------------------------------------------------------------
% Direct contrasts versus MCI
% -------------------------------------------------------------

Contrast = table(...
    ["CN vs MCI";"Dementia vs MCI"],...
    [C.Estimate(cn3);C.Estimate(dem3)],...
    [C.SE(cn3);C.SE(dem3)],...
    [C.Lower(cn3);C.Lower(dem3)],...
    [C.Upper(cn3);C.Upper(dem3)],...
    [C.pValue(cn3);C.pValue(dem3)],...
    'VariableNames',{...
    'Contrast','DifferenceInTauTime',...
    'SE','CI_Lower','CI_Upper','PValue'});

fprintf('\n============================================================\n');
fprintf(' DIFFERENCE IN TAU x TIME ASSOCIATION VS MCI\n');
fprintf('============================================================\n');

disp(Contrast);

%% ------------------------------------------------------------
% Save outputs
% -------------------------------------------------------------

% Information criteria are stored in ModelCriterion
ICred  = mdlReduced.ModelCriterion;
ICfull = mdlFull.ModelCriterion;

% Global interaction test table
GlobalTest = table(...
    LR,...
    dfDifference,...
    pGlobal,...
    ICred.AIC,...
    ICfull.AIC,...
    ICred.BIC,...
    ICfull.BIC,...
    'VariableNames',{...
    'LRStat',...
    'DF',...
    'PValue',...
    'ReducedAIC',...
    'FullAIC',...
    'ReducedBIC',...
    'FullBIC'});

% Save stage-specific estimates
writetable(StageSpecific,...
    fullfile(outputDir,...
    'Step42B_DiagnosisSpecific_TauTime_Estimates.csv'));

% Save pairwise contrasts
writetable(Contrast,...
    fullfile(outputDir,...
    'Step42B_TauTime_Diagnosis_Contrasts.csv'));

% Save global test
writetable(GlobalTest,...
    fullfile(outputDir,...
    'Step42B_Global_TauTimeDiagnosis_Test.csv'));

% Save complete models/workspace
save(fullfile(outputDir,...
    'Step42B_TauTimeDiagnosis_Models.mat'),...
    'mdlReduced','mdlFull','StageSpecific',...
    'Contrast','GlobalTest');

%% ------------------------------------------------------------
% Complete
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' STEP 42B COMPLETE\n');
fprintf('============================================================\n');

fprintf('Participants: 143: PASS\n');
fprintf('CN / MCI / Dementia: 33 / 75 / 35: PASS\n');
fprintf('Reference diagnosis: MCI\n');

fprintf('\nGLOBAL Tau PET x Time x Diagnosis:\n');
fprintf('LR chi-square(%d) = %.3f, P = %.6f\n',...
    dfDifference,LR,pGlobal);

fprintf('\nMODEL FIT\n');
fprintf('Reduced model: AIC = %.3f | BIC = %.3f\n',...
    ICred.AIC,ICred.BIC);
fprintf('Full model:    AIC = %.3f | BIC = %.3f\n',...
    ICfull.AIC,ICfull.BIC);

fprintf('\nSaved outputs to:\n%s\n',outputDir);

fprintf('\nFiles created:\n');
fprintf('  Step42B_DiagnosisSpecific_TauTime_Estimates.csv\n');
fprintf('  Step42B_TauTime_Diagnosis_Contrasts.csv\n');
fprintf('  Step42B_Global_TauTimeDiagnosis_Test.csv\n');
fprintf('  Step42B_TauTimeDiagnosis_Models.mat\n');

fprintf('\n============================================================\n');