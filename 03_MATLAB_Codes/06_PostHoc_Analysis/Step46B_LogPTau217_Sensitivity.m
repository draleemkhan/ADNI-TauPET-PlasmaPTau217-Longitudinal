%% STEP 46B: LOG-TRANSFORMED p-tau217 SENSITIVITY
% Tests whether biomarker entry-order findings are sensitive to
% log transformation of plasma p-tau217.
%
% Sequence 1:
% Centiloid -> log p-tau217 -> tau PET
%
% Sequence 2:
% Centiloid -> tau PET -> log p-tau217

clear; clc;

outputDir = 'D:\output';

fprintf('\n============================================================\n');
fprintf(' STEP 46B: LOG p-tau217 SENSITIVITY ANALYSIS\n');
fprintf('============================================================\n');

%% ------------------------------------------------------------
% Load locked model-ready dataset
% -------------------------------------------------------------

T = readtable(fullfile(outputDir,...
    'Step08_ModelReady_379_Observations.csv'),...
    'VariableNamingRule','preserve');

fprintf('\nDataset: %d observations\n',height(T));
fprintf('Participants: %d\n',numel(unique(T.RID)));

assert(height(T)==379,'Expected 379 observations.');
assert(numel(unique(T.RID))==143,'Expected 143 participants.');

fprintf('Dataset audit: PASS\n');

%% ------------------------------------------------------------
% Check raw p-tau217
% -------------------------------------------------------------

assert(ismember('PTau217_Fujirebio',T.Properties.VariableNames),...
    'Raw p-tau217 variable missing.');

assert(all(~isnan(T.PTau217_Fujirebio)),...
    'Missing raw p-tau217 values detected.');

assert(all(T.PTau217_Fujirebio > 0),...
    'Non-positive p-tau217 value detected; log transform invalid.');

fprintf('\nRaw p-tau217 range: %.6f to %.6f\n',...
    min(T.PTau217_Fujirebio),...
    max(T.PTau217_Fujirebio));

%% ------------------------------------------------------------
% Natural-log transformation
% -------------------------------------------------------------

T.LogPTau217 = log(T.PTau217_Fujirebio);

%% ------------------------------------------------------------
% Standardise using one value per participant
%
% Biomarker is baseline/static, so participant-level
% standardisation avoids weighting participants according to
% number of longitudinal CDR-SB observations.
% -------------------------------------------------------------

[~,ia] = unique(T.RID,'stable');

participantLogPTau = T.LogPTau217(ia);

muLogPTau = mean(participantLogPTau,'omitnan');
sdLogPTau = std(participantLogPTau,'omitnan');

T.LogPTau217_z = (T.LogPTau217 - muLogPTau) ./ sdLogPTau;

fprintf('\nLOG p-tau217 STANDARDISATION\n');
fprintf('Participant-level N: %d\n',numel(participantLogPTau));
fprintf('Mean log p-tau217: %.6f\n',muLogPTau);
fprintf('SD log p-tau217:   %.6f\n',sdLogPTau);

fprintf('\nRaw p-tau217:\n');
fprintf('  Median = %.4f\n',median(T.PTau217_Fujirebio));
fprintf('  Min    = %.4f\n',min(T.PTau217_Fujirebio));
fprintf('  Max    = %.4f\n',max(T.PTau217_Fujirebio));

fprintf('\nLog p-tau217:\n');
fprintf('  Median = %.4f\n',median(T.LogPTau217));
fprintf('  Min    = %.4f\n',min(T.LogPTau217));
fprintf('  Max    = %.4f\n',max(T.LogPTau217));

%% ------------------------------------------------------------
% Sex coding identical to primary analysis
% -------------------------------------------------------------

T.SexModel = categorical(string(T.Sex),...
    {'Male','Female'});

fprintf('\nSex reference category: Male\n');

%% ------------------------------------------------------------
% Model formulas
% -------------------------------------------------------------

% Model 1: Centiloid
f1 = ['CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
      'Centiloid_z + TimeYears:Centiloid_z + ' ...
      '(1|RID) + (TimeYears-1|RID)'];

% Model 2L: Centiloid + log p-tau217
f2L = ['CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
       'Centiloid_z + LogPTau217_z + ' ...
       'TimeYears:Centiloid_z + TimeYears:LogPTau217_z + ' ...
       '(1|RID) + (TimeYears-1|RID)'];

% Model 2T: Centiloid + tau PET
f2T = ['CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
       'Centiloid_z + Tau_z + ' ...
       'TimeYears:Centiloid_z + TimeYears:Tau_z + ' ...
       '(1|RID) + (TimeYears-1|RID)'];

% Full model
f3 = ['CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
      'Centiloid_z + LogPTau217_z + Tau_z + ' ...
      'TimeYears:Centiloid_z + TimeYears:LogPTau217_z + ' ...
      'TimeYears:Tau_z + ' ...
      '(1|RID) + (TimeYears-1|RID)'];

%% ------------------------------------------------------------
% Fit ML models
% -------------------------------------------------------------

fprintf('\nFitting Model 1: Centiloid...\n');
mdl1 = fitlme(T,f1,...
    'FitMethod','ML',...
    'DummyVarCoding','reference');

fprintf('Fitting Model 2L: + log p-tau217...\n');
mdl2L = fitlme(T,f2L,...
    'FitMethod','ML',...
    'DummyVarCoding','reference');

fprintf('Fitting Model 2T: + tau PET...\n');
mdl2T = fitlme(T,f2T,...
    'FitMethod','ML',...
    'DummyVarCoding','reference');

fprintf('Fitting Model 3: full model...\n');
mdl3 = fitlme(T,f3,...
    'FitMethod','ML',...
    'DummyVarCoding','reference');

%% ------------------------------------------------------------
% Forward sequence
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' FORWARD: CENTILOID -> LOG p-tau217 -> TAU PET\n');
fprintf('============================================================\n');

disp(compare(mdl1,mdl2L));
disp(compare(mdl2L,mdl3));

%% ------------------------------------------------------------
% Reverse sequence
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' REVERSE: CENTILOID -> TAU PET -> LOG p-tau217\n');
fprintf('============================================================\n');

disp(compare(mdl1,mdl2T));
disp(compare(mdl2T,mdl3));

%% ------------------------------------------------------------
% Model criteria
% -------------------------------------------------------------

c1  = mdl1.ModelCriterion;
c2L = mdl2L.ModelCriterion;
c2T = mdl2T.ModelCriterion;
c3  = mdl3.ModelCriterion;

Model = [
    "Model 1"
    "Model 2L"
    "Model 2T"
    "Model 3"
    ];

BiomarkerBlock = [
    "Centiloid"
    "Centiloid + log plasma p-tau217"
    "Centiloid + tau PET"
    "Centiloid + log plasma p-tau217 + tau PET"
    ];

AIC = [c1.AIC;c2L.AIC;c2T.AIC;c3.AIC];
BIC = [c1.BIC;c2L.BIC;c2T.BIC;c3.BIC];

LogLikelihood = [
    mdl1.LogLikelihood
    mdl2L.LogLikelihood
    mdl2T.LogLikelihood
    mdl3.LogLikelihood
    ];

ModelFit = table(Model,BiomarkerBlock,AIC,BIC,LogLikelihood);

fprintf('\n============================================================\n');
fprintf(' MODEL FIT\n');
fprintf('============================================================\n');

disp(ModelFit);

%% ------------------------------------------------------------
% Incremental tests
% -------------------------------------------------------------

LR_L_after_C  = 2*(mdl2L.LogLikelihood-mdl1.LogLikelihood);
LR_T_after_CL = 2*(mdl3.LogLikelihood-mdl2L.LogLikelihood);

LR_T_after_C  = 2*(mdl2T.LogLikelihood-mdl1.LogLikelihood);
LR_L_after_CT = 2*(mdl3.LogLikelihood-mdl2T.LogLikelihood);

df = 2;

P_L_after_C  = 1-chi2cdf(LR_L_after_C,df);
P_T_after_CL = 1-chi2cdf(LR_T_after_CL,df);
P_T_after_C  = 1-chi2cdf(LR_T_after_C,df);
P_L_after_CT = 1-chi2cdf(LR_L_after_CT,df);

Addition = [
    "Log p-tau217 after Centiloid"
    "Tau PET after Centiloid + log p-tau217"
    "Tau PET after Centiloid"
    "Log p-tau217 after Centiloid + tau PET"
    ];

LRStat = [
    LR_L_after_C
    LR_T_after_CL
    LR_T_after_C
    LR_L_after_CT
    ];

DF = repmat(2,4,1);

PValue = [
    P_L_after_C
    P_T_after_CL
    P_T_after_C
    P_L_after_CT
    ];

DeltaAIC = [
    c2L.AIC-c1.AIC
    c3.AIC-c2L.AIC
    c2T.AIC-c1.AIC
    c3.AIC-c2T.AIC
    ];

DeltaBIC = [
    c2L.BIC-c1.BIC
    c3.BIC-c2L.BIC
    c2T.BIC-c1.BIC
    c3.BIC-c2T.BIC
    ];

IncrementalTests = table(...
    Addition,LRStat,DF,PValue,DeltaAIC,DeltaBIC);

fprintf('\n============================================================\n');
fprintf(' INCREMENTAL TESTS\n');
fprintf('============================================================\n');

disp(IncrementalTests);

%% ------------------------------------------------------------
% Full model fixed effects
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' FULL LOG-p-tau217 MODEL FIXED EFFECTS\n');
fprintf('============================================================\n');

disp(mdl3.Coefficients);

%% ------------------------------------------------------------
% Sequence consistency
% -------------------------------------------------------------

forwardTotal = LR_L_after_C + LR_T_after_CL;
reverseTotal = LR_T_after_C + LR_L_after_CT;

fprintf('\n============================================================\n');
fprintf(' SEQUENCE CONSISTENCY CHECK\n');
fprintf('============================================================\n');

fprintf('Forward total LR: %.8f\n',forwardTotal);
fprintf('Reverse total LR: %.8f\n',reverseTotal);
fprintf('Absolute difference: %.12f\n',...
    abs(forwardTotal-reverseTotal));

assert(abs(forwardTotal-reverseTotal)<1e-6,...
    'Sequence consistency failure.');

fprintf('Sequence consistency: PASS\n');

%% ------------------------------------------------------------
% Save
% -------------------------------------------------------------

writetable(ModelFit,...
    fullfile(outputDir,...
    'Step46B_LogPTau217_ModelComparison.csv'));

writetable(IncrementalTests,...
    fullfile(outputDir,...
    'Step46B_LogPTau217_IncrementalTests.csv'));

TransformSummary = table(...
    muLogPTau,sdLogPTau,...
    min(T.PTau217_Fujirebio),...
    median(T.PTau217_Fujirebio),...
    max(T.PTau217_Fujirebio),...
    'VariableNames',{...
    'MeanLogPTau217','SDLogPTau217',...
    'RawMin','RawMedian','RawMax'});

writetable(TransformSummary,...
    fullfile(outputDir,...
    'Step46B_LogPTau217_TransformationSummary.csv'));

save(fullfile(outputDir,...
    'Step46B_LogPTau217_Models.mat'),...
    'mdl1','mdl2L','mdl2T','mdl3',...
    'ModelFit','IncrementalTests',...
    'muLogPTau','sdLogPTau');

%% ------------------------------------------------------------
% Final summary
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' STEP 46B COMPLETE\n');
fprintf('============================================================\n');

fprintf('Observations: 379: PASS\n');
fprintf('Participants: 143: PASS\n');

fprintf('\nLOG p-tau217 AFTER CENTILOID\n');
fprintf('LR chi-square(2) = %.4f | P = %.6g\n',...
    LR_L_after_C,P_L_after_C);
fprintf('Delta AIC = %.4f | Delta BIC = %.4f\n',...
    c2L.AIC-c1.AIC,c2L.BIC-c1.BIC);

fprintf('\nTAU PET AFTER CENTILOID + LOG p-tau217\n');
fprintf('LR chi-square(2) = %.4f | P = %.6g\n',...
    LR_T_after_CL,P_T_after_CL);
fprintf('Delta AIC = %.4f | Delta BIC = %.4f\n',...
    c3.AIC-c2L.AIC,c3.BIC-c2L.BIC);

fprintf('\nTAU PET AFTER CENTILOID\n');
fprintf('LR chi-square(2) = %.4f | P = %.6g\n',...
    LR_T_after_C,P_T_after_C);
fprintf('Delta AIC = %.4f | Delta BIC = %.4f\n',...
    c2T.AIC-c1.AIC,c2T.BIC-c1.BIC);

fprintf('\nLOG p-tau217 AFTER CENTILOID + TAU PET\n');
fprintf('LR chi-square(2) = %.4f | P = %.6g\n',...
    LR_L_after_CT,P_L_after_CT);
fprintf('Delta AIC = %.4f | Delta BIC = %.4f\n',...
    c3.AIC-c2T.AIC,c3.BIC-c2T.BIC);

fprintf('\nSaved outputs to:\n%s\n',outputDir);
fprintf('============================================================\n');