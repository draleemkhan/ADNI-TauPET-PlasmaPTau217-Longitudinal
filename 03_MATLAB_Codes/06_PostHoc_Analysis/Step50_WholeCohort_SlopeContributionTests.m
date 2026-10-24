%% ============================================================
% STEP 50
% WHOLE-COHORT SINGLE AND JOINT TAU-BIOMARKER SLOPE TESTS
%
% Purpose:
% Whole-cohort analogue of validated Step 49I MCI analysis.
%
% Tests:
% 1. Tau PET x Time beyond Centiloid
% 2. p-tau217 x Time beyond Centiloid
% 3. Joint Tau PET + p-tau217 slopes beyond Centiloid
%
% IMPORTANT:
% - Whole analytic cohort: expected 143 participants / 379 observations
% - ML used for nested fixed-effect comparisons
% - Random effects:
%       uncorrelated random intercept + random TimeYears slope
%       (1|RID) + (-1 + TimeYears|RID)
%
% - Single-biomarker tests deliberately exclude the OTHER tau biomarker.
% - Joint test retains BOTH tau-biomarker main effects and removes
%   only their two Time interaction terms.
%% ============================================================

clear;
clc;

fprintf('\n============================================================\n');
fprintf(' STEP 50: WHOLE-COHORT SLOPE CONTRIBUTION TESTS\n');
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
% Sex variable exactly as Step 47
% -------------------------------------------------------------

if ismember('Sex',T.Properties.VariableNames)
    sexSource = 'Sex';
elseif ismember('Sex_cat',T.Properties.VariableNames)
    sexSource = 'Sex_cat';
else
    error('Sex variable not found.');
end

T.SexModel = categorical(string(T.(sexSource)));

T.SexModel = reordercats(T.SexModel,...
    {'Male','Female'});

fprintf('Sex source: %s\n',sexSource);
fprintf('Sex reference category: Male\n');

%% ------------------------------------------------------------
% Whole-cohort audit
% -------------------------------------------------------------

nObs = height(T);
nParticipants = numel(unique(T.RID));

fprintf('\nWhole-cohort observations: %d\n',nObs);
fprintf('Whole-cohort participants: %d\n',nParticipants);

assert(nObs==379,...
    'Expected 379 observations.');

assert(nParticipants==143,...
    'Expected 143 participants.');

fprintf('Whole-cohort audit: PASS\n');

%% ============================================================
% STEP 47 WHOLE-COHORT REPRODUCTION GATE
%% ============================================================

fprintf('\n============================================================\n');
fprintf(' STEP 47 WHOLE-COHORT REPRODUCTION GATE\n');
fprintf('============================================================\n');

S47 = load(step47File,'mdlFull');

assert(isfield(S47,'mdlFull'),...
    'mdlFull not found in Step47 MAT file.');

savedModel = S47.mdlFull;

savedCrit = savedModel.ModelCriterion;

savedAIC = savedCrit.AIC;
savedBIC = savedCrit.BIC;
savedLL  = savedModel.LogLikelihood;

savedCoef = savedModel.Coefficients;
savedNames = string(savedCoef.Name);

tauTerm  = 'TimeYears:Tau_z';
ptauTerm = 'TimeYears:PTau217_z';

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

%% ------------------------------------------------------------
% Explicit reproduction formula
% -------------------------------------------------------------

fullFormula = [...
    'CDRSB ~ 1 + Age + Education + SexModel + '...
    'TimeYears*Tau_z + '...
    'TimeYears*Centiloid_z + '...
    'TimeYears*PTau217_z + '...
    '(1|RID) + (-1 + TimeYears|RID)'];

fprintf('\nCandidate reproduction formula:\n%s\n\n',fullFormula);

fprintf('Fitting candidate whole-cohort reproduction model using ML...\n');

mdlRepro = fitlme(T,fullFormula,...
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
% Reproduction comparison
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
    fprintf(' STEP 47 WHOLE-COHORT REPRODUCTION: FAIL\n');
    fprintf('============================================================\n');

    error(['Step47 whole-cohort model was not exactly reproduced. '...
           'Slope tests aborted.']);
end

fprintf('\n============================================================\n');
fprintf(' STEP 47 WHOLE-COHORT REPRODUCTION: PASS\n');
fprintf('============================================================\n');

%% ============================================================
% TEST 1
% TAU PET x TIME BEYOND CENTILOID
%
% p-tau217 deliberately excluded.
% Tau and Centiloid main effects retained.
%% ============================================================

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

fprintf('\nFitting Tau PET-only whole-cohort models...\n');

mdlTauReduced = fitlme(T,tauReducedFormula,...
    'FitMethod','ML');

mdlTauFull = fitlme(T,tauFullFormula,...
    'FitMethod','ML');

cmpTau = compare(mdlTauReduced,mdlTauFull);

%% ============================================================
% TEST 2
% p-tau217 x TIME BEYOND CENTILOID
%
% Tau PET deliberately excluded.
% p-tau217 and Centiloid main effects retained.
%% ============================================================

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

fprintf('Fitting p-tau217-only whole-cohort models...\n');

mdlPTauReduced = fitlme(T,ptauReducedFormula,...
    'FitMethod','ML');

mdlPTauFull = fitlme(T,ptauFullFormula,...
    'FitMethod','ML');

cmpPTau = compare(mdlPTauReduced,mdlPTauFull);

%% ============================================================
% TEST 3
% JOINT TAU PET + p-tau217 SLOPES BEYOND CENTILOID
%
% Both tau biomarker main effects retained.
% Reduced model removes ONLY:
%   Tau PET x Time
%   p-tau217 x Time
%
% Centiloid x Time retained.
%% ============================================================

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

fprintf('Fitting joint tau-biomarker whole-cohort models...\n');

mdlJointReduced = fitlme(T,jointReducedFormula,...
    'FitMethod','ML');

mdlJointFull = fitlme(T,jointFullFormula,...
    'FitMethod','ML');

cmpJoint = compare(mdlJointReduced,mdlJointFull);

%% ============================================================
% DISPLAY LIKELIHOOD-RATIO TESTS
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
% Extract LRT statistics
% -------------------------------------------------------------

[tauLR,tauDF,tauP,tauDAIC,tauDBIC] = ...
    getComparisonStats(mdlTauReduced,mdlTauFull);

[ptauLR,ptauDF,ptauP,ptauDAIC,ptauDBIC] = ...
    getComparisonStats(mdlPTauReduced,mdlPTauFull);

[jointLR,jointDF,jointP,jointDAIC,jointDBIC] = ...
    getComparisonStats(mdlJointReduced,mdlJointFull);

%% ------------------------------------------------------------
% Single-biomarker slope coefficients
% -------------------------------------------------------------

coefTau = mdlTauFull.Coefficients;
namesTau = string(coefTau.Name);

idxTau = find(namesTau==tauTerm,1);

assert(~isempty(idxTau),...
    'Tau PET x Time coefficient not found.');

tauBeta = coefTau.Estimate(idxTau);
tauSE = coefTau.SE(idxTau);
tauLower = coefTau.Lower(idxTau);
tauUpper = coefTau.Upper(idxTau);
tauWaldP = coefTau.pValue(idxTau);

coefPTau = mdlPTauFull.Coefficients;
namesPTau = string(coefPTau.Name);

idxPTau = find(namesPTau==ptauTerm,1);

assert(~isempty(idxPTau),...
    'p-tau217 x Time coefficient not found.');

ptauBeta = coefPTau.Estimate(idxPTau);
ptauSE = coefPTau.SE(idxPTau);
ptauLower = coefPTau.Lower(idxPTau);
ptauUpper = coefPTau.Upper(idxPTau);
ptauWaldP = coefPTau.pValue(idxPTau);

fprintf('\n============================================================\n');
fprintf(' SINGLE-BIOMARKER WHOLE-COHORT SLOPE ESTIMATES\n');
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
% Final result table
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
fprintf(' STEP 50 FINAL RESULTS\n');
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
    'Step50_WholeCohort_SlopeContributionTests.csv'));

writetable(SingleBiomarkerSlopes,...
    fullfile(outputDir,...
    'Step50_WholeCohort_SingleBiomarker_SlopeCoefficients.csv'));

save(fullfile(outputDir,...
    'Step50_WholeCohort_SlopeContributionModels.mat'),...
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
fprintf(' STEP 50 COMPLETE\n');
fprintf('============================================================\n');

fprintf('Whole-cohort participants: %d: PASS\n',nParticipants);
fprintf('Whole-cohort observations: %d: PASS\n',nObs);
fprintf('Step 47 whole-cohort reproduction: PASS\n');

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
fprintf('  Step50_WholeCohort_SlopeContributionTests.csv\n');
fprintf('  Step50_WholeCohort_SingleBiomarker_SlopeCoefficients.csv\n');
fprintf('  Step50_WholeCohort_SlopeContributionModels.mat\n');

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

% These comparisons differ only in fixed effects.
% Derive df from the number of fixed-effect coefficients.
df = full.NumCoefficients - reduced.NumCoefficients;

assert(df > 0,...
    'Invalid likelihood-ratio test degrees of freedom.');

p = 1-chi2cdf(LR,df);

dAIC = full.ModelCriterion.AIC - ...
    reduced.ModelCriterion.AIC;

dBIC = full.ModelCriterion.BIC - ...
    reduced.ModelCriterion.BIC;

end