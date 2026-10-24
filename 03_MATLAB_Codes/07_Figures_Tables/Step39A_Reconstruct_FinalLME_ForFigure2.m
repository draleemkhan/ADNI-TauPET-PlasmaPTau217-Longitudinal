%% Step39A_Reconstruct_FinalLME_ForFigure2.m
% ============================================================
% FIGURE 2 — RECONSTRUCT AND VERIFY FINAL PRIMARY LME
% Corrected sex coding:
%   Source = Sex
%   Reference = Male
% ============================================================

clear; clc;

fprintf('\n============================================================\n');
fprintf(' STEP 39A: RECONSTRUCT FINAL PRIMARY LME FOR FIGURE 2\n');
fprintf('============================================================\n');

outputDir = ...
    'D:\output';

%% ============================================================
% LOAD MODEL-READY DATA
% ============================================================

dataFile = fullfile(outputDir,...
    'Step08_ModelReady_379_Observations.csv');

T = readtable(dataFile,...
    'VariableNamingRule','preserve',...
    'TextType','string');

fprintf('\nDataset: %d observations, %d participants\n',...
    height(T),numel(unique(T.RID)));

assert(height(T)==379,...
    'Expected 379 observations.');

assert(numel(unique(T.RID))==143,...
    'Expected 143 participants.');

%% ============================================================
% REQUIRED VARIABLES
% ============================================================

required = [
    "RID"
    "CDRSB"
    "TimeYears"
    "Age"
    "Education"
    "Tau_z"
    "Centiloid_z"
    "PTau217_z"
    "Sex"
];

vn = string(T.Properties.VariableNames);

for k=1:numel(required)

    assert(ismember(required(k),vn),...
        'Missing required variable: %s',required(k));

end

fprintf('Required-variable audit: PASS\n');

%% ============================================================
% DATA TYPES
% ============================================================

T.RID = categorical(T.RID);

% Explicitly reconstruct the exact model sex factor.
T.SexModel = categorical(T.Sex);

sexCats = categories(T.SexModel);

assert(ismember('Male',sexCats),...
    'Male category not found.');

assert(ismember('Female',sexCats),...
    'Female category not found.');

% Force Male to reference level.
T.SexModel = reordercats(T.SexModel,...
    {'Male','Female'});

fprintf('\nSex coding used for reconstructed model:\n');
fprintf('Reference: Male\n');
fprintf('Comparison: Female vs Male\n');

disp(categories(T.SexModel));

%% ============================================================
% ENSURE NUMERIC VARIABLES ARE NUMERIC
% ============================================================

numericVars = [
    "CDRSB"
    "TimeYears"
    "Age"
    "Education"
    "Tau_z"
    "Centiloid_z"
    "PTau217_z"
];

for k=1:numel(numericVars)

    v = numericVars(k);

    if ~isnumeric(T.(v))
        T.(v) = str2double(T.(v));
    end

    assert(all(isfinite(T.(v))),...
        'Nonfinite values detected in %s.',v);

end

%% ============================================================
% LOCKED FINAL PRIMARY MODEL
%
% Random intercept + random TimeYears slope
% with ZERO intercept-slope covariance.
%
% Separate random-effect terms enforce the uncorrelated structure.
% ============================================================

formula = [ ...
    'CDRSB ~ TimeYears + Age + Education + ' ...
    'Tau_z + Centiloid_z + PTau217_z + SexModel + ' ...
    'TimeYears:Tau_z + ' ...
    'TimeYears:Centiloid_z + ' ...
    'TimeYears:PTau217_z + ' ...
    '(1|RID) + (-1 + TimeYears|RID)' ...
    ];

fprintf('\nFitting locked primary model using REML...\n');

lme = fitlme(T,formula,...
    'FitMethod','REML');

fprintf('\nRECONSTRUCTED FIXED EFFECTS\n');
disp(lme.Coefficients);

%% ============================================================
% AUTHORITATIVE STEP 10 RESULTS
% ============================================================

authFile = fullfile(outputDir,...
    'Step10_FinalPrimary_UncorrelatedLME_FixedEffects.csv');

A = readtable(authFile,...
    'VariableNamingRule','preserve',...
    'TextType','string');

% Convert authoritative numeric fields if necessary.
numAuth = ["Estimate","SE","Lower","Upper","pValue"];

for k=1:numel(numAuth)
    if ~isnumeric(A.(numAuth(k)))
        A.(numAuth(k)) = str2double(A.(numAuth(k)));
    end
end

fprintf('\nAUTHORITATIVE STEP 10 FIXED EFFECTS\n');
disp(A);

%% ============================================================
% NORMALISE TERM NAMES
% ============================================================

newNames = string(lme.CoefficientNames)';
authNames = string(A.Name);

newNames  = replace(newNames,"\:",":");
authNames = replace(authNames,"\:",":");

% MATLAB should generate SexModel_Female because Male is reference.
fprintf('\nRECONSTRUCTED TERM NAMES\n');
disp(newNames);

assert(any(newNames=="SexModel_Female"),...
    ['Expected SexModel_Female term was not generated. ' ...
     'Check sex reference coding.']);

%% ============================================================
% MATCH EVERY FIXED EFFECT
% ============================================================

assert(numel(newNames)==height(A),...
    'Number of fixed effects differs from authoritative Step 10.');

comparison = table;

comparison.Term = authNames;
comparison.AuthoritativeBeta = A.Estimate;
comparison.ReconstructedBeta = nan(height(A),1);
comparison.AbsoluteDifference = nan(height(A),1);

for k=1:height(A)

    idx = find(newNames==authNames(k));

    assert(numel(idx)==1,...
        'Could not uniquely match term: %s',authNames(k));

    comparison.ReconstructedBeta(k) = ...
        lme.Coefficients.Estimate(idx);

    comparison.AbsoluteDifference(k) = ...
        abs(comparison.AuthoritativeBeta(k) - ...
        comparison.ReconstructedBeta(k));

end

fprintf('\n============================================================\n');
fprintf(' COEFFICIENT REPRODUCTION CHECK\n');
fprintf('============================================================\n');

disp(comparison);

maxDiff = max(comparison.AbsoluteDifference);

fprintf('Maximum absolute coefficient difference: %.12g\n',...
    maxDiff);

%% ============================================================
% VERIFY REPRODUCTION
% ============================================================

tol = 1e-6;

if maxDiff <= tol

    fprintf('\nFIXED-EFFECT REPRODUCTION: PASS\n');

else

    fprintf('\nFIXED-EFFECT REPRODUCTION: FAIL\n');
    fprintf(['Do NOT generate Figure 2 yet. ' ...
             'We will reconcile the model specification first.\n']);
    return;

end

%% ============================================================
% TAU × TIME BENCHMARK
% ============================================================

tauIdx = find(authNames=="TimeYears:Tau_z");

assert(numel(tauIdx)==1,...
    'Could not uniquely identify Tau PET x Time.');

tauBeta = comparison.ReconstructedBeta(tauIdx);

fprintf('\nTau PET x Time reconstructed beta: %.6f\n',...
    tauBeta);

assert(abs(tauBeta-0.413141)<0.001,...
    'Tau PET x Time benchmark failed.');

fprintf('Tau PET x Time benchmark: PASS\n');

%% ============================================================
% FIXED-EFFECT COVARIANCE MATRIX
% ============================================================

CovBeta = lme.CoefficientCovariance;
coefNames = string(lme.CoefficientNames)';

assert(size(CovBeta,1)==numel(coefNames));
assert(size(CovBeta,2)==numel(coefNames));

fprintf('\nFixed-effect covariance matrix: %d x %d\n',...
    size(CovBeta,1),size(CovBeta,2));

%% ============================================================
% SAVE
% ============================================================

saveFile = fullfile(outputDir,...
    'Step39A_FinalPrimaryLME_ForFigure2.mat');

save(saveFile,...
    'lme',...
    'CovBeta',...
    'coefNames',...
    'comparison',...
    '-v7.3');

covTable = array2table(CovBeta,...
    'VariableNames',...
    matlab.lang.makeValidName(cellstr(coefNames)),...
    'RowNames',...
    cellstr(coefNames));

writetable(covTable,...
    fullfile(outputDir,...
    'Step39A_FinalPrimary_FixedEffectCovariance.csv'),...
    'WriteRowNames',true);

writetable(comparison,...
    fullfile(outputDir,...
    'Step39A_FinalPrimary_CoefficientVerification.csv'));

%% ============================================================
% COMPLETE
% ============================================================

fprintf('\n============================================================\n');
fprintf(' STEP 39A COMPLETE\n');
fprintf('============================================================\n');

fprintf('Dataset: 379 observations / 143 participants: PASS\n');
fprintf('Sex reference = Male: PASS\n');
fprintf('Fixed-effect reproduction: PASS\n');
fprintf('Tau PET x Time benchmark: PASS\n');

fprintf('\nSaved model/covariance source:\n%s\n',saveFile);

fprintf('\nREADY FOR FIGURE 2.\n');