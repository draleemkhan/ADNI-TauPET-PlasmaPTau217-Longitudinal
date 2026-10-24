%% ============================================================
% STEP 44: PLASMA p-tau217 x TIME x DIAGNOSIS INTERACTION
%
% Purpose:
%   Formally test whether the association between baseline
%   plasma p-tau217 and longitudinal CDR-SB change differs
%   according to baseline diagnostic stage.
%
% Mirrors Step 42B for Tau PET.
%
% Cohort:
%   143 participants
%   379 longitudinal observations
%
% Diagnosis:
%   Cognitively normal = 33
%   MCI                = 75
%   Dementia           = 35
%
% Reference diagnosis:
%   MCI
%
% Global test:
%   Reduced model versus full model using ML.
%   Full model adds:
%       PTau217_z x TimeYears x Diagnosis
%
% Outputs:
%   Step44_DiagnosisSpecific_PTauTime_Estimates.csv
%   Step44_PTauTime_Diagnosis_Contrasts.csv
%   Step44_Global_PTauTimeDiagnosis_Test.csv
%   Step44_PTauTimeDiagnosis_Models.mat
% ============================================================

clear;
clc;

fprintf('\n============================================================\n');
fprintf(' STEP 44: p-tau217 x TIME x DIAGNOSIS INTERACTION\n');
fprintf('============================================================\n');

%% ------------------------------------------------------------
% Paths
% -------------------------------------------------------------

outputDir = ...
    'D:\output';

dataFile = fullfile(outputDir,...
    'Step08_ModelReady_379_Observations.csv');

%% ------------------------------------------------------------
% Load data
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
    "Sex"
    "Baseline_Diagnosis"];

for k = 1:numel(required)

    assert(ismember(required(k),...
        string(T.Properties.VariableNames)),...
        'Missing required variable: %s',required(k));

end

nParticipants = numel(unique(T.RID));

assert(height(T)==379,...
    'Expected 379 observations.');

assert(nParticipants==143,...
    'Expected 143 participants.');

fprintf('Participants: %d\n',nParticipants);
fprintf('Dataset audit: PASS\n');

%% ------------------------------------------------------------
% Sex coding
% Male reference
% -------------------------------------------------------------

T.SexModel = categorical(string(T.Sex));

assert(any(strcmp(categories(T.SexModel),'Male')),...
    'Male sex category not found.');

assert(any(strcmp(categories(T.SexModel),'Female')),...
    'Female sex category not found.');

T.SexModel = reordercats(T.SexModel,...
    {'Male','Female'});

%% ------------------------------------------------------------
% Diagnosis coding
% MCI reference
% -------------------------------------------------------------

T.DxModel = categorical(...
    string(T.Baseline_Diagnosis));

dxCats = categories(T.DxModel);

assert(any(strcmp(dxCats,'MCI')),...
    'MCI diagnosis category not found.');

assert(any(strcmp(dxCats,'Cognitively normal')),...
    'Cognitively normal category not found.');

assert(any(strcmp(dxCats,'Dementia')),...
    'Dementia category not found.');

T.DxModel = reordercats(T.DxModel,...
    {'MCI','Cognitively normal','Dementia'});

fprintf('\nDiagnosis reference category: MCI\n');

%% ------------------------------------------------------------
% Participant-level diagnosis audit
% -------------------------------------------------------------

[~,ia] = unique(T.RID,'stable');

P = T(ia,{'RID','Baseline_Diagnosis'});

dxParticipant = string(P.Baseline_Diagnosis);

nMCI = sum(dxParticipant=="MCI");
nCN  = sum(dxParticipant=="Cognitively normal");
nDEM = sum(dxParticipant=="Dementia");

fprintf('\nParticipant-level diagnosis counts:\n');
fprintf('MCI:                 %d\n',nMCI);
fprintf('Cognitively normal:  %d\n',nCN);
fprintf('Dementia:            %d\n',nDEM);

assert(nMCI==75,...
    'Expected 75 MCI participants.');

assert(nCN==33,...
    'Expected 33 cognitively normal participants.');

assert(nDEM==35,...
    'Expected 35 dementia participants.');

fprintf('Diagnosis-count audit: PASS\n');

%% ------------------------------------------------------------
% MODEL SPECIFICATION
%
% Reduced:
%   Diagnosis main effects
%   Diagnosis x time
%
%   Tau PET main + Tau PET x time
%   Centiloid main + Centiloid x time
%   p-tau217 main + p-tau217 x time
%
%   No diagnosis-specific modification of p-tau217 x time.
%
% Full:
%   Adds:
%     p-tau217 x diagnosis
%     p-tau217 x time x diagnosis
%
% The global LRT therefore tests whether p-tau217 x time
% differs across diagnostic stages while preserving hierarchy.
% -------------------------------------------------------------

formulaReduced = [...
    'CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
    'Tau_z + Centiloid_z + PTau217_z + DxModel + ' ...
    'TimeYears:Tau_z + ' ...
    'TimeYears:Centiloid_z + ' ...
    'TimeYears:PTau217_z + ' ...
    'TimeYears:DxModel + ' ...
    'PTau217_z:DxModel + ' ...
    '(1|RID) + (TimeYears-1|RID)'];

formulaFull = [...
    'CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
    'Tau_z + Centiloid_z + PTau217_z + DxModel + ' ...
    'TimeYears:Tau_z + ' ...
    'TimeYears:Centiloid_z + ' ...
    'TimeYears:PTau217_z + ' ...
    'TimeYears:DxModel + ' ...
    'PTau217_z:DxModel + ' ...
    'TimeYears:PTau217_z:DxModel + ' ...
    '(1|RID) + (TimeYears-1|RID)'];

%% ------------------------------------------------------------
% Fit ML models
% -------------------------------------------------------------

fprintf('\nFitting reduced ML model...\n');

mdlReduced = fitlme(T,formulaReduced,...
    'FitMethod','ML');

fprintf('\nFitting full ML model...\n');

mdlFull = fitlme(T,formulaFull,...
    'FitMethod','ML');

%% ------------------------------------------------------------
% Display full-model coefficients
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' FULL MODEL FIXED EFFECTS\n');
fprintf('============================================================\n\n');

disp(mdlFull.Coefficients);

%% ------------------------------------------------------------
% Global likelihood-ratio test
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' GLOBAL p-tau217 x TIME x DIAGNOSIS TEST\n');
fprintf('============================================================\n\n');

cmp = compare(mdlReduced,mdlFull);

disp(cmp);

LLred  = mdlReduced.LogLikelihood;
LLfull = mdlFull.LogLikelihood;

LR = 2*(LLfull-LLred);

dfDifference = ...
    mdlFull.NumCoefficients - mdlReduced.NumCoefficients;

pGlobal = 1-chi2cdf(LR,dfDifference);

fprintf('\nReduced log-likelihood: %.6f\n',LLred);
fprintf('Full log-likelihood:    %.6f\n',LLfull);
fprintf('LR statistic:           %.6f\n',LR);
fprintf('df difference:          %d\n',dfDifference);
fprintf('GLOBAL interaction P:   %.8g\n',pGlobal);

%% ------------------------------------------------------------
% Identify required coefficient names
% -------------------------------------------------------------

C = mdlFull.Coefficients;

coefNames = string(C.Name);

baseName = "TimeYears:PTau217_z";

cn3Name = ...
    "TimeYears:PTau217_z:DxModel_Cognitively normal";

dem3Name = ...
    "TimeYears:PTau217_z:DxModel_Dementia";

idxBase = find(coefNames==baseName);
idxCN3  = find(coefNames==cn3Name);
idxDEM3 = find(coefNames==dem3Name);

% MATLAB may print interactions in another equivalent ordering.
% If exact strings are not found, identify terms by components.

if isempty(idxBase)

    idxBase = find(...
        contains(coefNames,"TimeYears") & ...
        contains(coefNames,"PTau217_z") & ...
        ~contains(coefNames,"DxModel"));

end

if isempty(idxCN3)

    idxCN3 = find(...
        contains(coefNames,"TimeYears") & ...
        contains(coefNames,"PTau217_z") & ...
        contains(coefNames,"DxModel_Cognitively normal"));

end

if isempty(idxDEM3)

    idxDEM3 = find(...
        contains(coefNames,"TimeYears") & ...
        contains(coefNames,"PTau217_z") & ...
        contains(coefNames,"DxModel_Dementia"));

end

assert(numel(idxBase)==1,...
    'Could not uniquely identify baseline p-tau217 x time term.');

assert(numel(idxCN3)==1,...
    'Could not uniquely identify CN three-way interaction.');

assert(numel(idxDEM3)==1,...
    'Could not uniquely identify dementia three-way interaction.');

fprintf('\nIdentified p-tau217/time interaction terms:\n');
disp(coefNames([idxBase idxCN3 idxDEM3]));

%% ------------------------------------------------------------
% Fixed-effect covariance matrix
% -------------------------------------------------------------

V = mdlFull.CoefficientCovariance;

b = C.Estimate;

%% ------------------------------------------------------------
% Diagnosis-specific p-tau217 x time estimates
%
% MCI:
%   beta(Time x p-tau217)
%
% CN:
%   beta(Time x p-tau217)
%   + beta(Time x p-tau217 x CN)
%
% Dementia:
%   beta(Time x p-tau217)
%   + beta(Time x p-tau217 x Dementia)
% -------------------------------------------------------------

p = height(C);

L_MCI = zeros(1,p);
L_CN  = zeros(1,p);
L_DEM = zeros(1,p);

L_MCI(idxBase) = 1;

L_CN(idxBase) = 1;
L_CN(idxCN3)  = 1;

L_DEM(idxBase) = 1;
L_DEM(idxDEM3) = 1;

% Estimates
betaMCI = L_MCI*b;
betaCN  = L_CN*b;
betaDEM = L_DEM*b;

% SE
seMCI = sqrt(L_MCI*V*L_MCI');
seCN  = sqrt(L_CN*V*L_CN');
seDEM = sqrt(L_DEM*V*L_DEM');

% Residual DF
dfResid = mdlFull.DFE;

crit = tinv(0.975,dfResid);

% CI
loMCI = betaMCI - crit*seMCI;
hiMCI = betaMCI + crit*seMCI;

loCN = betaCN - crit*seCN;
hiCN = betaCN + crit*seCN;

loDEM = betaDEM - crit*seDEM;
hiDEM = betaDEM + crit*seDEM;

% P values
pMCI = 2*(1-tcdf(abs(betaMCI/seMCI),dfResid));
pCN  = 2*(1-tcdf(abs(betaCN/seCN),dfResid));
pDEM = 2*(1-tcdf(abs(betaDEM/seDEM),dfResid));

%% ------------------------------------------------------------
% Diagnosis-specific table
% -------------------------------------------------------------

Diagnosis = [...
    "MCI"
    "Cognitively normal"
    "Dementia"];

N = [...
    nMCI
    nCN
    nDEM];

Beta = [...
    betaMCI
    betaCN
    betaDEM];

SE = [...
    seMCI
    seCN
    seDEM];

CI_Lower = [...
    loMCI
    loCN
    loDEM];

CI_Upper = [...
    hiMCI
    hiCN
    hiDEM];

PValue = [...
    pMCI
    pCN
    pDEM];

StageSpecific = table(...
    Diagnosis,...
    N,...
    Beta,...
    SE,...
    CI_Lower,...
    CI_Upper,...
    PValue);

fprintf('\n============================================================\n');
fprintf(' DIAGNOSIS-SPECIFIC p-tau217 x TIME ESTIMATES\n');
fprintf(' FROM FORMAL INTERACTION MODEL\n');
fprintf('============================================================\n\n');

disp(StageSpecific);

%% ------------------------------------------------------------
% Contrasts versus MCI
%
% CN vs MCI = three-way CN coefficient
% Dementia vs MCI = three-way dementia coefficient
% -------------------------------------------------------------

L_CNvsMCI = zeros(1,p);
L_CNvsMCI(idxCN3) = 1;

L_DEMvsMCI = zeros(1,p);
L_DEMvsMCI(idxDEM3) = 1;

diffCN = L_CNvsMCI*b;
seDiffCN = sqrt(L_CNvsMCI*V*L_CNvsMCI');

diffDEM = L_DEMvsMCI*b;
seDiffDEM = sqrt(L_DEMvsMCI*V*L_DEMvsMCI');

loDiffCN = diffCN - crit*seDiffCN;
hiDiffCN = diffCN + crit*seDiffCN;

loDiffDEM = diffDEM - crit*seDiffDEM;
hiDiffDEM = diffDEM + crit*seDiffDEM;

pDiffCN = 2*(1-tcdf(abs(diffCN/seDiffCN),dfResid));
pDiffDEM = 2*(1-tcdf(abs(diffDEM/seDiffDEM),dfResid));

ContrastName = [...
    "CN vs MCI"
    "Dementia vs MCI"];

DifferenceInPTauTime = [...
    diffCN
    diffDEM];

ContrastSE = [...
    seDiffCN
    seDiffDEM];

ContrastCI_Lower = [...
    loDiffCN
    loDiffDEM];

ContrastCI_Upper = [...
    hiDiffCN
    hiDiffDEM];

ContrastPValue = [...
    pDiffCN
    pDiffDEM];

Contrast = table(...
    ContrastName,...
    DifferenceInPTauTime,...
    ContrastSE,...
    ContrastCI_Lower,...
    ContrastCI_Upper,...
    ContrastPValue);

fprintf('\n============================================================\n');
fprintf(' DIFFERENCE IN p-tau217 x TIME ASSOCIATION VS MCI\n');
fprintf('============================================================\n\n');

disp(Contrast);

%% ------------------------------------------------------------
% Save model-fit information
% -------------------------------------------------------------

ICred  = mdlReduced.ModelCriterion;
ICfull = mdlFull.ModelCriterion;

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

%% ------------------------------------------------------------
% Save outputs
% -------------------------------------------------------------

writetable(StageSpecific,...
    fullfile(outputDir,...
    'Step44_DiagnosisSpecific_PTauTime_Estimates.csv'));

writetable(Contrast,...
    fullfile(outputDir,...
    'Step44_PTauTime_Diagnosis_Contrasts.csv'));

writetable(GlobalTest,...
    fullfile(outputDir,...
    'Step44_Global_PTauTimeDiagnosis_Test.csv'));

save(fullfile(outputDir,...
    'Step44_PTauTimeDiagnosis_Models.mat'),...
    'mdlReduced','mdlFull',...
    'StageSpecific','Contrast',...
    'GlobalTest');

%% ------------------------------------------------------------
% Complete
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' STEP 44 COMPLETE\n');
fprintf('============================================================\n');

fprintf('Participants: 143: PASS\n');
fprintf('CN / MCI / Dementia: 33 / 75 / 35: PASS\n');
fprintf('Reference diagnosis: MCI\n');

fprintf('\nGLOBAL p-tau217 x Time x Diagnosis:\n');
fprintf('LR chi-square(%d) = %.3f, P = %.6f\n',...
    dfDifference,LR,pGlobal);

fprintf('\nDIAGNOSIS-SPECIFIC ESTIMATES\n');
fprintf('CN:       beta = %.3f, 95%% CI %.3f to %.3f, P = %.4f\n',...
    betaCN,loCN,hiCN,pCN);

fprintf('MCI:      beta = %.3f, 95%% CI %.3f to %.3f, P = %.4f\n',...
    betaMCI,loMCI,hiMCI,pMCI);

fprintf('Dementia: beta = %.3f, 95%% CI %.3f to %.3f, P = %.4f\n',...
    betaDEM,loDEM,hiDEM,pDEM);

fprintf('\nMODEL FIT\n');
fprintf('Reduced: AIC = %.3f | BIC = %.3f\n',...
    ICred.AIC,ICred.BIC);

fprintf('Full:    AIC = %.3f | BIC = %.3f\n',...
    ICfull.AIC,ICfull.BIC);

fprintf('\nSaved outputs:\n');
fprintf('  Step44_DiagnosisSpecific_PTauTime_Estimates.csv\n');
fprintf('  Step44_PTauTime_Diagnosis_Contrasts.csv\n');
fprintf('  Step44_Global_PTauTimeDiagnosis_Test.csv\n');
fprintf('  Step44_PTauTimeDiagnosis_Models.mat\n');

fprintf('\n============================================================\n');