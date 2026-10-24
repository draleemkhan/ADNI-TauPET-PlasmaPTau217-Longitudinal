%% ============================================================
% STEP 49I
% EXACT RANDOM-EFFECTS REPRODUCTION + MCI SLOPE TESTS
%
% Purpose:
% 1. Reproduce the SAVED Step 47 MCI model using an explicit
%    uncorrelated random-intercept/random-slope specification.
%
% 2. Only if Step 47 reproduction PASSES, test:
%    A. Tau PET x Time beyond Centiloid
%    B. p-tau217 x Time beyond Centiloid
%    C. Joint Tau PET + p-tau217 slopes beyond Centiloid
%
% IMPORTANT:
% - ML is used for nested fixed-effect comparisons.
% - The reproduction gate is mandatory.
% - No inferential results should be interpreted if the gate fails.
%% ============================================================

clear;
clc;

fprintf('\n============================================================\n');
fprintf(' STEP 49I: EXACT RANDOM-EFFECTS MCI SLOPE TESTS\n');
fprintf('============================================================\n');

%% ------------------------------------------------------------
% Paths
% -------------------------------------------------------------

outputDir = ...
    'D:\output';

dataFile = fullfile(outputDir,...
    'Step08_ModelReady_379_Observations.csv');

step47File = fullfile(outputDir,...
    'Step47_LevelVsSlope_Models.mat');

assert(isfile(dataFile),...
    'Step08 model-ready dataset not found.');

assert(isfile(step47File),...
    'Step47 model file not found.');

%% ------------------------------------------------------------
% Load data
% -------------------------------------------------------------

T = readtable(dataFile,...
    'VariableNamingRule','preserve');

fprintf('\nFull dataset: %d observations\n',height(T));

required = {...
    'RID',...
    'CDRSB',...
    'TimeYears',...
    'Age',...
    'Education',...
    'Tau_z',...
    'Centiloid_z',...
    'PTau217_z'};

for k = 1:numel(required)
    assert(ismember(required{k},T.Properties.VariableNames),...
        'Missing required variable: %s',required{k});
end

fprintf('Required-variable audit: PASS\n');

%% ------------------------------------------------------------
% Identify sex and diagnosis variables exactly as Step 47
% -------------------------------------------------------------

if ismember('Sex',T.Properties.VariableNames)
    sexSource = 'Sex';
elseif ismember('Sex_cat',T.Properties.VariableNames)
    sexSource = 'Sex_cat';
else
    error('Sex variable not found.');
end

if ismember('Baseline_Diagnosis',T.Properties.VariableNames)
    dxSource = 'Baseline_Diagnosis';
elseif ismember('Diagnosis',T.Properties.VariableNames)
    dxSource = 'Diagnosis';
else
    error('Diagnosis variable not found.');
end

T.SexModel = categorical(string(T.(sexSource)));
T.SexModel = reordercats(T.SexModel,...
    {'Male','Female'});

T.DxModel = categorical(string(T.(dxSource)));

fprintf('Sex source: %s\n',sexSource);
fprintf('Diagnosis source: %s\n',dxSource);
fprintf('Sex reference category: Male\n');

%% ------------------------------------------------------------
% MCI subset
% -------------------------------------------------------------

isMCI = string(T.(dxSource))=="MCI";
TMCI = T(isMCI,:);

nObsMCI = height(TMCI);
nMCI = numel(unique(TMCI.RID));

fprintf('\nMCI observations: %d\n',nObsMCI);
fprintf('MCI participants: %d\n',nMCI);

assert(nObsMCI==208,...
    'Expected 208 MCI observations.');

assert(nMCI==75,...
    'Expected 75 MCI participants.');

fprintf('MCI cohort audit: PASS\n');

%% ------------------------------------------------------------
% Load saved Step 47 model
% -------------------------------------------------------------

S47 = load(step47File,'mdlFullMCI');

assert(isfield(S47,'mdlFullMCI'),...
    'mdlFullMCI not found in Step47 MAT file.');

savedModel = S47.mdlFullMCI;

savedCrit = savedModel.ModelCriterion;

savedAIC = savedCrit.AIC;
savedBIC = savedCrit.BIC;
savedLL  = savedModel.LogLikelihood;

savedCoef = savedModel.Coefficients;

tauTerm = 'TimeYears:Tau_z';
ptauTerm = 'TimeYears:PTau217_z';

savedNames = string(savedCoef.Name);

idxTauSaved = find(savedNames==tauTerm,1);
idxPTauSaved = find(savedNames==ptauTerm,1);

assert(~isempty(idxTauSaved),...
    'Tau PET x Time term not found in saved Step47 model.');

assert(~isempty(idxPTauSaved),...
    'p-tau217 x Time term not found in saved Step47 model.');

savedTauBeta = savedCoef.Estimate(idxTauSaved);
savedTauSE   = savedCoef.SE(idxTauSaved);

savedPTauBeta = savedCoef.Estimate(idxPTauSaved);
savedPTauSE   = savedCoef.SE(idxPTauSaved);

%% ============================================================
% REPRODUCTION MODEL
%% ============================================================

fprintf('\n============================================================\n');
fprintf(' STEP 47 REPRODUCTION GATE\n');
fprintf('============================================================\n');

% Exact fixed effects from Step 47.
%
% Critical change:
% Random intercept and random TimeYears slope are specified as
% separate terms so that no intercept-slope covariance is fitted.

fullFormula = [...
    'CDRSB ~ 1 + Age + Education + SexModel + '...
    'TimeYears*Tau_z + '...
    'TimeYears*Centiloid_z + '...
    'TimeYears*PTau217_z + '...
    '(1|RID) + (-1 + TimeYears|RID)'];

fprintf('\nCandidate reproduction formula:\n%s\n\n',fullFormula);

fprintf('Fitting candidate reproduction model using ML...\n');

mdlRepro = fitlme(TMCI,fullFormula,...
    'FitMethod','ML');

newCrit = mdlRepro.ModelCriterion;

newAIC = newCrit.AIC;
newBIC = newCrit.BIC;
newLL  = mdlRepro.LogLikelihood;

newCoef = mdlRepro.Coefficients;
newNames = string(newCoef.Name);

idxTauNew = find(newNames==tauTerm,1);
idxPTauNew = find(newNames==ptauTerm,1);

assert(~isempty(idxTauNew),...
    'Tau PET x Time term not found in reproduction model.');

assert(~isempty(idxPTauNew),...
    'p-tau217 x Time term not found in reproduction model.');

newTauBeta = newCoef.Estimate(idxTauNew);
newTauSE   = newCoef.SE(idxTauNew);

newPTauBeta = newCoef.Estimate(idxPTauNew);
newPTauSE   = newCoef.SE(idxPTauNew);

%% ------------------------------------------------------------
% Display reproduction comparison
% -------------------------------------------------------------

fprintf('\nSaved Step47 AIC: %.12f\n',savedAIC);
fprintf('New model AIC:    %.12f\n',newAIC);

fprintf('\nSaved Step47 BIC: %.12f\n',savedBIC);
fprintf('New model BIC:    %.12f\n',newBIC);

fprintf('\nSaved Step47 LogLik: %.12f\n',savedLL);
fprintf('New model LogLik:    %.12f\n',newLL);

fprintf('\nTau PET x Time\n');
fprintf('Saved beta = %.9f | New beta = %.9f\n',...
    savedTauBeta,newTauBeta);
fprintf('Saved SE   = %.9f | New SE   = %.9f\n',...
    savedTauSE,newTauSE);

fprintf('\np-tau217 x Time\n');
fprintf('Saved beta = %.9f | New beta = %.9f\n',...
    savedPTauBeta,newPTauBeta);
fprintf('Saved SE   = %.9f | New SE   = %.9f\n',...
    savedPTauSE,newPTauSE);

%% ------------------------------------------------------------
% Numerical reproduction gate
% -------------------------------------------------------------

tol = 1e-6;

passAIC = abs(newAIC-savedAIC) < tol;
passBIC = abs(newBIC-savedBIC) < tol;
passLL  = abs(newLL-savedLL) < tol;

passTauBeta = abs(newTauBeta-savedTauBeta) < tol;
passTauSE   = abs(newTauSE-savedTauSE) < tol;

passPTauBeta = abs(newPTauBeta-savedPTauBeta) < tol;
passPTauSE   = abs(newPTauSE-savedPTauSE) < tol;

fprintf('\nREPRODUCTION CHECK\n');

fprintf('AIC:                %s\n',passFail(passAIC));
fprintf('BIC:                %s\n',passFail(passBIC));
fprintf('LogLikelihood:      %s\n',passFail(passLL));
fprintf('Tau beta:           %s\n',passFail(passTauBeta));
fprintf('Tau SE:             %s\n',passFail(passTauSE));
fprintf('p-tau217 beta:      %s\n',passFail(passPTauBeta));
fprintf('p-tau217 SE:        %s\n',passFail(passPTauSE));

reproductionPASS = ...
    passAIC && passBIC && passLL && ...
    passTauBeta && passTauSE && ...
    passPTauBeta && passPTauSE;

if ~reproductionPASS

    fprintf('\n============================================================\n');
    fprintf(' STEP 47 REPRODUCTION: FAIL\n');
    fprintf('============================================================\n');

    error(['Explicit uncorrelated random-effects model did not '...
        'reproduce Step47. Definitive slope tests aborted.']);
end

fprintf('\n============================================================\n');
fprintf(' STEP 47 REPRODUCTION: PASS\n');
fprintf('============================================================\n');

%% ============================================================
% DEFINITIVE MCI SLOPE TESTS
%% ============================================================

fprintf('\nProceeding to definitive MCI slope tests...\n');

%% ------------------------------------------------------------
% TEST 1
% Tau PET x Time beyond Centiloid
%
% p-tau217 is NOT included.
% Tau and Centiloid main effects retained.
% Reduced model removes ONLY Tau x Time.
% -------------------------------------------------------------

tauReducedFormula = [...
    'CDRSB ~ 1 + Age + Education + SexModel + '...
    'TimeYears + Tau_z + Centiloid_z + '...
    'TimeYears:Centiloid_z + '...
    '(1|RID) + (-1 + TimeYears|RID)'];

tauFullFormula = [...
    'CDRSB ~ 1 + Age + Education + SexModel + '...
    'TimeYears + Tau_z + Centiloid_z + '...
    'TimeYears:Tau_z + TimeYears:Centiloid_z + '...
    '(1|RID) + (-1 + TimeYears|RID)'];

fprintf('\nFitting Tau PET-only slope models...\n');

mdlTauReduced = fitlme(TMCI,tauReducedFormula,...
    'FitMethod','ML');

mdlTauFull = fitlme(TMCI,tauFullFormula,...
    'FitMethod','ML');

cmpTau = compare(mdlTauReduced,mdlTauFull);

%% ------------------------------------------------------------
% TEST 2
% p-tau217 x Time beyond Centiloid
%
% Tau PET is NOT included.
% p-tau217 and Centiloid main effects retained.
% Reduced model removes ONLY p-tau217 x Time.
% -------------------------------------------------------------

ptauReducedFormula = [...
    'CDRSB ~ 1 + Age + Education + SexModel + '...
    'TimeYears + PTau217_z + Centiloid_z + '...
    'TimeYears:Centiloid_z + '...
    '(1|RID) + (-1 + TimeYears|RID)'];

ptauFullFormula = [...
    'CDRSB ~ 1 + Age + Education + SexModel + '...
    'TimeYears + PTau217_z + Centiloid_z + '...
    'TimeYears:PTau217_z + TimeYears:Centiloid_z + '...
    '(1|RID) + (-1 + TimeYears|RID)'];

fprintf('Fitting p-tau217-only slope models...\n');

mdlPTauReduced = fitlme(TMCI,ptauReducedFormula,...
    'FitMethod','ML');

mdlPTauFull = fitlme(TMCI,ptauFullFormula,...
    'FitMethod','ML');

cmpPTau = compare(mdlPTauReduced,mdlPTauFull);

%% ------------------------------------------------------------
% TEST 3
% Joint Tau PET + p-tau217 slopes beyond Centiloid
%
% ALL biomarker main effects retained in BOTH models.
% Reduced model removes ONLY:
%   Tau PET x Time
%   p-tau217 x Time
%
% Centiloid x Time retained.
% -------------------------------------------------------------

jointReducedFormula = [...
    'CDRSB ~ 1 + Age + Education + SexModel + '...
    'TimeYears + Tau_z + Centiloid_z + PTau217_z + '...
    'TimeYears:Centiloid_z + '...
    '(1|RID) + (-1 + TimeYears|RID)'];

jointFullFormula = [...
    'CDRSB ~ 1 + Age + Education + SexModel + '...
    'TimeYears + Tau_z + Centiloid_z + PTau217_z + '...
    'TimeYears:Tau_z + '...
    'TimeYears:Centiloid_z + '...
    'TimeYears:PTau217_z + '...
    '(1|RID) + (-1 + TimeYears|RID)'];

fprintf('Fitting joint tau-biomarker slope models...\n');

mdlJointReduced = fitlme(TMCI,jointReducedFormula,...
    'FitMethod','ML');

mdlJointFull = fitlme(TMCI,jointFullFormula,...
    'FitMethod','ML');

cmpJoint = compare(mdlJointReduced,mdlJointFull);

%% ============================================================
% DISPLAY LRTs
%% ============================================================

fprintf('\n============================================================\n');
fprintf(' TEST 1: TAU PET x TIME BEYOND CENTILOID\n');
fprintf('============================================================\n');
disp(cmpTau);

fprintf('\n============================================================\n');
fprintf(' TEST 2: p-tau217 x TIME BEYOND CENTILOID\n');
fprintf('============================================================\n');
disp(cmpPTau);

fprintf('\n============================================================\n');
fprintf(' TEST 3: JOINT TAU-BIOMARKER SLOPES BEYOND CENTILOID\n');
fprintf('============================================================\n');
disp(cmpJoint);

%% ------------------------------------------------------------
% Extract LRT statistics directly from likelihoods
% -------------------------------------------------------------

[tauLR,tauDF,tauP,tauDAIC,tauDBIC] = ...
    getComparisonStats(mdlTauReduced,mdlTauFull);

[ptauLR,ptauDF,ptauP,ptauDAIC,ptauDBIC] = ...
    getComparisonStats(mdlPTauReduced,mdlPTauFull);

[jointLR,jointDF,jointP,jointDAIC,jointDBIC] = ...
    getComparisonStats(mdlJointReduced,mdlJointFull);

%% ------------------------------------------------------------
% Extract single-biomarker slope coefficients
% -------------------------------------------------------------

coefTau = mdlTauFull.Coefficients;
namesTau = string(coefTau.Name);

idx = find(namesTau==tauTerm,1);

tauBeta = coefTau.Estimate(idx);
tauSE = coefTau.SE(idx);
tauLower = coefTau.Lower(idx);
tauUpper = coefTau.Upper(idx);
tauWaldP = coefTau.pValue(idx);

coefPTau = mdlPTauFull.Coefficients;
namesPTau = string(coefPTau.Name);

idx = find(namesPTau==ptauTerm,1);

ptauBeta = coefPTau.Estimate(idx);
ptauSE = coefPTau.SE(idx);
ptauLower = coefPTau.Lower(idx);
ptauUpper = coefPTau.Upper(idx);
ptauWaldP = coefPTau.pValue(idx);

fprintf('\n============================================================\n');
fprintf(' SINGLE-BIOMARKER MCI SLOPE ESTIMATES\n');
fprintf('============================================================\n');

fprintf('\nTau PET x Time beyond Centiloid:\n');
fprintf('beta = %.6f\n',tauBeta);
fprintf('SE = %.6f\n',tauSE);
fprintf('95%% CI = %.6f to %.6f\n',tauLower,tauUpper);
fprintf('Wald P = %.8g\n',tauWaldP);

fprintf('\np-tau217 x Time beyond Centiloid:\n');
fprintf('beta = %.6f\n',ptauBeta);
fprintf('SE = %.6f\n',ptauSE);
fprintf('95%% CI = %.6f to %.6f\n',ptauLower,ptauUpper);
fprintf('Wald P = %.8g\n',ptauWaldP);

%% ------------------------------------------------------------
% Final results table
% -------------------------------------------------------------

TestName = [...
    "Tau PET x Time beyond Centiloid";...
    "p-tau217 x Time beyond Centiloid";...
    "Joint Tau PET + p-tau217 slopes beyond Centiloid"];

LRStat = [tauLR;ptauLR;jointLR];
DF = [tauDF;ptauDF;jointDF];
PValue = [tauP;ptauP;jointP];
DeltaAIC = [tauDAIC;ptauDAIC;jointDAIC];
DeltaBIC = [tauDBIC;ptauDBIC;jointDBIC];

FinalResults = table(...
    TestName,LRStat,DF,PValue,DeltaAIC,DeltaBIC);

fprintf('\n============================================================\n');
fprintf(' STEP 49I FINAL RESULTS\n');
fprintf('============================================================\n');

disp(FinalResults);

%% ------------------------------------------------------------
% Coefficient table
% -------------------------------------------------------------

Biomarker = ["Tau PET";"Plasma p-tau217"];

Beta = [tauBeta;ptauBeta];
SE = [tauSE;ptauSE];
CI_Lower = [tauLower;ptauLower];
CI_Upper = [tauUpper;ptauUpper];
WaldPValue = [tauWaldP;ptauWaldP];

SingleBiomarkerSlopes = table(...
    Biomarker,Beta,SE,CI_Lower,CI_Upper,WaldPValue);

%% ------------------------------------------------------------
% Save outputs
% -------------------------------------------------------------

writetable(FinalResults,...
    fullfile(outputDir,...
    'Step49I_MCI_SlopeContributionTests.csv'));

writetable(SingleBiomarkerSlopes,...
    fullfile(outputDir,...
    'Step49I_MCI_SingleBiomarker_SlopeCoefficients.csv'));

save(fullfile(outputDir,...
    'Step49I_MCI_SlopeContributionModels.mat'),...
    'mdlRepro',...
    'mdlTauReduced','mdlTauFull',...
    'mdlPTauReduced','mdlPTauFull',...
    'mdlJointReduced','mdlJointFull',...
    'FinalResults','SingleBiomarkerSlopes',...
    'reproductionPASS');

%% ------------------------------------------------------------
% Final console summary
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' STEP 49I COMPLETE\n');
fprintf('============================================================\n');

fprintf('MCI participants: %d: PASS\n',nMCI);
fprintf('MCI observations: %d: PASS\n',nObsMCI);
fprintf('Step 47 reproduction: PASS\n');

fprintf('\nPRIMARY SLOPE QUESTIONS\n');

fprintf('\n1. Tau PET x Time beyond Centiloid\n');
fprintf('LR chi-square(%d) = %.4f, P = %.8g\n',...
    tauDF,tauLR,tauP);
fprintf('Delta AIC = %.4f | Delta BIC = %.4f\n',...
    tauDAIC,tauDBIC);

fprintf('\n2. p-tau217 x Time beyond Centiloid\n');
fprintf('LR chi-square(%d) = %.4f, P = %.8g\n',...
    ptauDF,ptauLR,ptauP);
fprintf('Delta AIC = %.4f | Delta BIC = %.4f\n',...
    ptauDAIC,ptauDBIC);

fprintf('\n3. Joint Tau PET + p-tau217 slopes beyond Centiloid\n');
fprintf('LR chi-square(%d) = %.4f, P = %.8g\n',...
    jointDF,jointLR,jointP);
fprintf('Delta AIC = %.4f | Delta BIC = %.4f\n',...
    jointDAIC,jointDBIC);

fprintf('\nSaved outputs to:\n%s\n',outputDir);

fprintf('\nFiles created:\n');
fprintf('  Step49I_MCI_SlopeContributionTests.csv\n');
fprintf('  Step49I_MCI_SingleBiomarker_SlopeCoefficients.csv\n');
fprintf('  Step49I_MCI_SlopeContributionModels.mat\n');

fprintf('============================================================\n');

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function txt = passFail(tf)

if tf
    txt = 'PASS';
else
    txt = 'FAIL';
end

end

function [LR,df,p,dAIC,dBIC] = getComparisonStats(reduced,full)

LR = 2*(full.LogLikelihood-reduced.LogLikelihood);

df = full.NumCoefficients + ...
    length(full.covarianceParameters) - ...
    reduced.NumCoefficients - ...
    length(reduced.covarianceParameters);

% More robust df fallback based on model criterion if needed.
if df <= 0
    df = round(...
        (full.ModelCriterion.AIC - reduced.ModelCriterion.AIC + LR)/2);
end

p = 1-chi2cdf(LR,df);

dAIC = full.ModelCriterion.AIC - ...
    reduced.ModelCriterion.AIC;

dBIC = full.ModelCriterion.BIC - ...
    reduced.ModelCriterion.BIC;

end