%% Step19_Marginal_Conditional_R2.m
% ADNI Tau PET revised cohort analysis
%
% STEP 19 — Marginal and conditional R-squared
%
% Fits three nested longitudinal LMMs using the locked:
%   N = 143 participants
%   379 CDR-SB observations
%   uncorrelated random intercept + random Time slope
%
% Model 1: Centiloid + demographics
% Model 2: Model 1 + plasma p-tau217
% Model 3: Model 2 + tau PET
%
% Because the scientific outcome is longitudinal progression, each biomarker
% included in a model is represented by BOTH its main effect and interaction
% with Time.
%
% Model 1:
% CDRSB ~ Time + Centiloid + Age + Sex + Education + Centiloid*Time
%
% Model 2:
% Model 1 + p-tau217 + p-tau217*Time
%
% Model 3:
% Model 2 + Tau + Tau*Time
%
% Nakagawa-style variance decomposition:
%   R2 marginal    = Var(fixed predictions) / total variance
%   R2 conditional = [Var(fixed predictions) + mean random-effect variance]
%                    / total variance
%
% For random intercept + random Time slope with diagonal covariance,
% observation-specific random variance is:
%   Var(b0) + Time^2 * Var(b1)
% and its mean across observations is used.
%
% Historical approximate benchmarks:
%   Centiloid + demographics: R2m ~0.272, R2c ~0.915
%   + p-tau217:              R2m ~0.362, R2c ~0.909
%   + tau PET:               R2m ~0.392, R2c ~0.909
%
% These are in-sample explained-variance measures, NOT predictive accuracy.

clear; clc; close all;

outputDir = 'D:\output';
dataFile = fullfile(outputDir,'Step08_ModelReady_379_Observations.csv');

assert(isfile(dataFile),'Missing Step08 model-ready dataset.');

T = readtable(dataFile,'VariableNamingRule','preserve');

fprintf('\n====================================================\n');
fprintf(' STEP 19 - MARGINAL AND CONDITIONAL R-SQUARED\n');
fprintf('====================================================\n');

%% Numeric conversion
numVars = {'RID','CDRSB','TimeYears','Tau_z','Centiloid_z','PTau217_z','Age','Education'};
for k = 1:numel(numVars)
    if ~isnumeric(T.(numVars{k}))
        T.(numVars{k}) = str2double(string(T.(numVars{k})));
    end
end

T.RIDcat = categorical(T.RID);
T.SexModel = categorical(string(T.Sex),["Male","Female"]);

complete = isfinite(T.CDRSB) & isfinite(T.TimeYears) & ...
           isfinite(T.Tau_z) & isfinite(T.Centiloid_z) & ...
           isfinite(T.PTau217_z) & isfinite(T.Age) & ...
           isfinite(T.Education) & ~isundefined(T.SexModel);

M = T(complete,:);

fprintf('Complete observations: %d [expected 379]\n',height(M));
fprintf('Complete participants: %d [expected 143]\n',numel(unique(M.RID)));

%% Nested formulas
f1 = ['CDRSB ~ TimeYears + Centiloid_z + Age + SexModel + Education + ' ...
      'Centiloid_z:TimeYears + (1 + TimeYears|RIDcat)'];

f2 = ['CDRSB ~ TimeYears + Centiloid_z + PTau217_z + Age + SexModel + Education + ' ...
      'Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
      '(1 + TimeYears|RIDcat)'];

f3 = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + Age + SexModel + Education + ' ...
      'Tau_z:TimeYears + Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
      '(1 + TimeYears|RIDcat)'];

labels = ["Centiloid + demographics"; ...
          "+ p-tau217"; ...
          "+ tau PET"];

formulas = {f1,f2,f3};
models = cell(3,1);

Result = table();

for k = 1:3
    fprintf('\n----------------------------------------------------\n');
    fprintf('MODEL %d: %s\n',k,labels(k));
    fprintf('----------------------------------------------------\n');

    mdl = fitlme(M,formulas{k}, ...
        'FitMethod','REML', ...
        'CovariancePattern','Diagonal');

    models{k} = mdl;

    % Fixed-only fitted values
    fixedPred = fitted(mdl,'Conditional',false);
    varFixed = var(fixedPred,1);

    % Extract random-effect covariance matrix robustly
    [psi,~,stats] = covarianceParameters(mdl);

    % For this model there is one grouping factor with diagonal
    % random intercept/random slope covariance.
    if iscell(psi)
        G = psi{1};
    else
        G = psi;
    end

    if size(G,1) < 2 || size(G,2) < 2
        error('Expected 2x2 random-effects covariance matrix.');
    end

    varRI = G(1,1);
    varRS = G(2,2);
    covIR = G(1,2);

    % Observation-specific random-effect variance:
    % z_i G z_i', z_i = [1 TimeYears]
    t = double(M.TimeYears);
    randVarObs = varRI + 2*t*covIR + (t.^2)*varRS;
    varRandom = mean(randVarObs);

    % Residual variance
    residualVar = mdl.MSE;

    totalVar = varFixed + varRandom + residualVar;

    R2m = varFixed / totalVar;
    R2c = (varFixed + varRandom) / totalVar;

    fprintf('Fixed-effect variance: %.6f\n',varFixed);
    fprintf('Mean random-effect variance: %.6f\n',varRandom);
    fprintf('Residual variance: %.6f\n',residualVar);
    fprintf('Marginal R-squared: %.6f\n',R2m);
    fprintf('Conditional R-squared: %.6f\n',R2c);

    row = table(labels(k),varFixed,varRandom,residualVar,R2m,R2c, ...
        mdl.ModelCriterion.AIC,mdl.ModelCriterion.BIC, ...
        'VariableNames',{'Model','FixedVariance','RandomVariance', ...
        'ResidualVariance','Marginal_R2','Conditional_R2','AIC','BIC'});

    Result = [Result;row]; %#ok<AGROW>
end

%% Incremental marginal R2
Result.DeltaMarginalR2 = [NaN; ...
    Result.Marginal_R2(2)-Result.Marginal_R2(1); ...
    Result.Marginal_R2(3)-Result.Marginal_R2(2)];

Result.DeltaMarginalR2_PercentagePoints = 100*Result.DeltaMarginalR2;

fprintf('\n====================================================\n');
fprintf(' R-SQUARED SUMMARY\n');
fprintf('====================================================\n');
disp(Result);

fprintf('\nINCREMENTAL FIXED-EFFECT EXPLAINED VARIANCE\n');
fprintf('Adding p-tau217: Delta marginal R2 = %.4f (%.2f percentage points)\n', ...
    Result.DeltaMarginalR2(2),Result.DeltaMarginalR2_PercentagePoints(2));

fprintf('Adding tau PET beyond Centiloid + p-tau217: Delta marginal R2 = %.4f (%.2f percentage points)\n', ...
    Result.DeltaMarginalR2(3),Result.DeltaMarginalR2_PercentagePoints(3));

fprintf('\nIMPORTANT: These are in-sample variance-explanation measures,\n');
fprintf('not out-of-sample predictive accuracy.\n');

%% Save
writetable(Result, ...
    fullfile(outputDir,'Step19_Marginal_Conditional_R2.csv'));

save(fullfile(outputDir,'Step19_Marginal_Conditional_R2_Models.mat'), ...
    'models','formulas','Result');

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');
