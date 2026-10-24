%% ============================================================
% STEP 47
% DECOMPOSITION OF BIOMARKER INCREMENTAL CONTRIBUTION
% BASELINE LEVEL VS LONGITUDINAL SLOPE
%
% Purpose:
% Separately test the 1-df contribution of:
%   1) biomarker main effect (baseline level component)
%   2) biomarker x time interaction (longitudinal slope component)
%
% Analyses:
%   A. Whole analytic cohort
%   B. MCI subgroup
%
% Biomarkers:
%   - Tau PET
%   - Plasma p-tau217
%
% The key slope tests retain ALL biomarker main effects and the
% other biomarker x time interactions, and remove only the target
% biomarker x time term.
%
% Models are fitted using maximum likelihood because nested
% fixed-effect models are compared by likelihood-ratio tests.
%% ============================================================

clear;
clc;

fprintf('\n============================================================\n');
fprintf(' STEP 47: LEVEL VS SLOPE DECOMPOSITION\n');
fprintf('============================================================\n');

%% ------------------------------------------------------------
% Paths
% -------------------------------------------------------------

outputDir = ...
    'D:\output';

dataFile = fullfile(outputDir,...
    'Step08_ModelReady_379_Observations.csv');

assert(isfile(dataFile),...
    'Model-ready dataset not found.');

%% ------------------------------------------------------------
% Load dataset
% -------------------------------------------------------------

T = readtable(dataFile,...
    'VariableNamingRule','preserve');

fprintf('\nFull dataset: %d observations\n',height(T));

%% ------------------------------------------------------------
% Required-variable audit
% -------------------------------------------------------------

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

% Sex variable
if ismember('Sex',T.Properties.VariableNames)
    sexSource = 'Sex';
elseif ismember('Sex_cat',T.Properties.VariableNames)
    sexSource = 'Sex_cat';
else
    error('Sex variable not found.');
end

% Diagnosis variable
if ismember('Baseline_Diagnosis',T.Properties.VariableNames)
    dxSource = 'Baseline_Diagnosis';
elseif ismember('Diagnosis',T.Properties.VariableNames)
    dxSource = 'Diagnosis';
else
    error('Diagnosis variable not found.');
end

%% ------------------------------------------------------------
% Create model variables
% -------------------------------------------------------------

T.SexModel = categorical(string(T.(sexSource)));

% Force Male reference
T.SexModel = reordercats(T.SexModel,...
    {'Male','Female'});

T.DxModel = categorical(string(T.(dxSource)));

fprintf('Required-variable audit: PASS\n');

%% ------------------------------------------------------------
% Participant audit
% -------------------------------------------------------------

nFull = numel(unique(T.RID));

fprintf('Full cohort participants: %d\n',nFull);

assert(height(T)==379,...
    'Expected 379 observations.');

assert(nFull==143,...
    'Expected 143 participants.');

fprintf('Full-cohort audit: PASS\n');

%% ------------------------------------------------------------
% MCI subset
% -------------------------------------------------------------

isMCI = string(T.(dxSource))=="MCI";

TMCI = T(isMCI,:);

nMCI = numel(unique(TMCI.RID));

fprintf('\nMCI observations: %d\n',height(TMCI));
fprintf('MCI participants: %d\n',nMCI);

assert(height(TMCI)==208,...
    'Expected 208 MCI observations.');

assert(nMCI==75,...
    'Expected 75 MCI participants.');

fprintf('MCI audit: PASS\n');

%% ============================================================
% MODEL FORMULAS
%
% FULL MODEL:
% demographics +
% Tau + Centiloid + p-tau217 +
% Time x Tau +
% Time x Centiloid +
% Time x p-tau217
%
% Random structure follows primary model:
% uncorrelated random intercept + random slope
%% ============================================================

fullFormula = [...
    'CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
    'Tau_z + Centiloid_z + PTau217_z + ' ...
    'TimeYears:Tau_z + ' ...
    'TimeYears:Centiloid_z + ' ...
    'TimeYears:PTau217_z + ' ...
    '(1|RID) + (-1 + TimeYears|RID)'];

%% ------------------------------------------------------------
% Reduced models for SLOPE tests
%
% Important:
% ALL biomarker main effects remain in the model.
% Only the target biomarker x time term is removed.
% -------------------------------------------------------------

noTauSlopeFormula = [...
    'CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
    'Tau_z + Centiloid_z + PTau217_z + ' ...
    'TimeYears:Centiloid_z + ' ...
    'TimeYears:PTau217_z + ' ...
    '(1|RID) + (-1 + TimeYears|RID)'];

noPTauSlopeFormula = [...
    'CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
    'Tau_z + Centiloid_z + PTau217_z + ' ...
    'TimeYears:Tau_z + ' ...
    'TimeYears:Centiloid_z + ' ...
    '(1|RID) + (-1 + TimeYears|RID)'];

%% ------------------------------------------------------------
% Reduced models for MAIN-EFFECT tests
%
% Hierarchical caveat:
% These models retain the corresponding interaction while
% removing its main effect. They are included only as a formal
% decomposition of the fitted biomarker block.
% -------------------------------------------------------------

noTauMainFormula = [...
    'CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
    'Centiloid_z + PTau217_z + ' ...
    'TimeYears:Tau_z + ' ...
    'TimeYears:Centiloid_z + ' ...
    'TimeYears:PTau217_z + ' ...
    '(1|RID) + (-1 + TimeYears|RID)'];

noPTauMainFormula = [...
    'CDRSB ~ TimeYears + Age + Education + SexModel + ' ...
    'Tau_z + Centiloid_z + ' ...
    'TimeYears:Tau_z + ' ...
    'TimeYears:Centiloid_z + ' ...
    'TimeYears:PTau217_z + ' ...
    '(1|RID) + (-1 + TimeYears|RID)'];

%% ============================================================
% WHOLE COHORT
%% ============================================================

fprintf('\n============================================================\n');
fprintf(' WHOLE COHORT\n');
fprintf('============================================================\n');

fprintf('\nFitting full ML model...\n');

mdlFull = fitlme(T,fullFormula,...
    'FitMethod','ML');

fprintf('Fitting reduced slope models...\n');

mdlNoTauSlope = fitlme(T,noTauSlopeFormula,...
    'FitMethod','ML');

mdlNoPTauSlope = fitlme(T,noPTauSlopeFormula,...
    'FitMethod','ML');

fprintf('Fitting reduced main-effect models...\n');

mdlNoTauMain = fitlme(T,noTauMainFormula,...
    'FitMethod','ML');

mdlNoPTauMain = fitlme(T,noPTauMainFormula,...
    'FitMethod','ML');

%% ------------------------------------------------------------
% Compare models
% -------------------------------------------------------------

cmpTauSlope = compare(mdlNoTauSlope,mdlFull);
cmpPTauSlope = compare(mdlNoPTauSlope,mdlFull);

cmpTauMain = compare(mdlNoTauMain,mdlFull);
cmpPTauMain = compare(mdlNoPTauMain,mdlFull);

fprintf('\nTAU PET SLOPE CONTRIBUTION\n');
disp(cmpTauSlope);

fprintf('\np-tau217 SLOPE CONTRIBUTION\n');
disp(cmpPTauSlope);

fprintf('\nTAU PET MAIN-EFFECT CONTRIBUTION\n');
disp(cmpTauMain);

fprintf('\np-tau217 MAIN-EFFECT CONTRIBUTION\n');
disp(cmpPTauMain);

%% ------------------------------------------------------------
% Extract whole-cohort statistics
% -------------------------------------------------------------

[LR_TauSlope,df_TauSlope,p_TauSlope,...
 dAIC_TauSlope,dBIC_TauSlope] = ...
    extractComparison(cmpTauSlope);

[LR_PTauSlope,df_PTauSlope,p_PTauSlope,...
 dAIC_PTauSlope,dBIC_PTauSlope] = ...
    extractComparison(cmpPTauSlope);

[LR_TauMain,df_TauMain,p_TauMain,...
 dAIC_TauMain,dBIC_TauMain] = ...
    extractComparison(cmpTauMain);

[LR_PTauMain,df_PTauMain,p_PTauMain,...
 dAIC_PTauMain,dBIC_PTauMain] = ...
    extractComparison(cmpPTauMain);

Whole = table(...
    ["Tau PET";"Tau PET";"Plasma p-tau217";"Plasma p-tau217"],...
    ["Main effect";"Time interaction";...
     "Main effect";"Time interaction"],...
    [LR_TauMain;LR_TauSlope;LR_PTauMain;LR_PTauSlope],...
    [df_TauMain;df_TauSlope;df_PTauMain;df_PTauSlope],...
    [p_TauMain;p_TauSlope;p_PTauMain;p_PTauSlope],...
    [dAIC_TauMain;dAIC_TauSlope;dAIC_PTauMain;dAIC_PTauSlope],...
    [dBIC_TauMain;dBIC_TauSlope;dBIC_PTauMain;dBIC_PTauSlope],...
    'VariableNames',{...
    'Biomarker','Component','LRStat','DF','PValue',...
    'DeltaAIC','DeltaBIC'});

fprintf('\nWHOLE-COHORT DECOMPOSITION\n');
disp(Whole);

%% ============================================================
% MCI SUBGROUP
%% ============================================================

fprintf('\n============================================================\n');
fprintf(' MCI SUBGROUP\n');
fprintf('============================================================\n');

fprintf('\nFitting full MCI ML model...\n');

mdlFullMCI = fitlme(TMCI,fullFormula,...
    'FitMethod','ML');

fprintf('Fitting reduced MCI slope models...\n');

mdlNoTauSlopeMCI = fitlme(TMCI,noTauSlopeFormula,...
    'FitMethod','ML');

mdlNoPTauSlopeMCI = fitlme(TMCI,noPTauSlopeFormula,...
    'FitMethod','ML');

fprintf('Fitting reduced MCI main-effect models...\n');

mdlNoTauMainMCI = fitlme(TMCI,noTauMainFormula,...
    'FitMethod','ML');

mdlNoPTauMainMCI = fitlme(TMCI,noPTauMainFormula,...
    'FitMethod','ML');

%% ------------------------------------------------------------
% MCI comparisons
% -------------------------------------------------------------

cmpTauSlopeMCI = compare(mdlNoTauSlopeMCI,mdlFullMCI);
cmpPTauSlopeMCI = compare(mdlNoPTauSlopeMCI,mdlFullMCI);

cmpTauMainMCI = compare(mdlNoTauMainMCI,mdlFullMCI);
cmpPTauMainMCI = compare(mdlNoPTauMainMCI,mdlFullMCI);

fprintf('\nMCI: TAU PET SLOPE CONTRIBUTION\n');
disp(cmpTauSlopeMCI);

fprintf('\nMCI: p-tau217 SLOPE CONTRIBUTION\n');
disp(cmpPTauSlopeMCI);

fprintf('\nMCI: TAU PET MAIN-EFFECT CONTRIBUTION\n');
disp(cmpTauMainMCI);

fprintf('\nMCI: p-tau217 MAIN-EFFECT CONTRIBUTION\n');
disp(cmpPTauMainMCI);

%% ------------------------------------------------------------
% Extract MCI statistics
% -------------------------------------------------------------

[LR_TauSlopeMCI,df_TauSlopeMCI,p_TauSlopeMCI,...
 dAIC_TauSlopeMCI,dBIC_TauSlopeMCI] = ...
    extractComparison(cmpTauSlopeMCI);

[LR_PTauSlopeMCI,df_PTauSlopeMCI,p_PTauSlopeMCI,...
 dAIC_PTauSlopeMCI,dBIC_PTauSlopeMCI] = ...
    extractComparison(cmpPTauSlopeMCI);

[LR_TauMainMCI,df_TauMainMCI,p_TauMainMCI,...
 dAIC_TauMainMCI,dBIC_TauMainMCI] = ...
    extractComparison(cmpTauMainMCI);

[LR_PTauMainMCI,df_PTauMainMCI,p_PTauMainMCI,...
 dAIC_PTauMainMCI,dBIC_PTauMainMCI] = ...
    extractComparison(cmpPTauMainMCI);

MCI = table(...
    ["Tau PET";"Tau PET";"Plasma p-tau217";"Plasma p-tau217"],...
    ["Main effect";"Time interaction";...
     "Main effect";"Time interaction"],...
    [LR_TauMainMCI;LR_TauSlopeMCI;...
     LR_PTauMainMCI;LR_PTauSlopeMCI],...
    [df_TauMainMCI;df_TauSlopeMCI;...
     df_PTauMainMCI;df_PTauSlopeMCI],...
    [p_TauMainMCI;p_TauSlopeMCI;...
     p_PTauMainMCI;p_PTauSlopeMCI],...
    [dAIC_TauMainMCI;dAIC_TauSlopeMCI;...
     dAIC_PTauMainMCI;dAIC_PTauSlopeMCI],...
    [dBIC_TauMainMCI;dBIC_TauSlopeMCI;...
     dBIC_PTauMainMCI;dBIC_PTauSlopeMCI],...
    'VariableNames',{...
    'Biomarker','Component','LRStat','DF','PValue',...
    'DeltaAIC','DeltaBIC'});

fprintf('\nMCI DECOMPOSITION\n');
disp(MCI);

%% ============================================================
% Extract coefficients for the two slope terms
%% ============================================================

coefFull = mdlFull.Coefficients;
coefMCI  = mdlFullMCI.Coefficients;

WholeSlopeCoefficients = getSlopeCoefficients(coefFull);
MCISlopeCoefficients   = getSlopeCoefficients(coefMCI);

fprintf('\n============================================================\n');
fprintf(' WHOLE-COHORT SLOPE COEFFICIENTS\n');
fprintf('============================================================\n');
disp(WholeSlopeCoefficients);

fprintf('\n============================================================\n');
fprintf(' MCI SLOPE COEFFICIENTS\n');
fprintf('============================================================\n');
disp(MCISlopeCoefficients);

%% ============================================================
% Combined output
%% ============================================================

Whole.Population = repmat("Whole cohort",height(Whole),1);
MCI.Population   = repmat("MCI",height(MCI),1);

FinalResults = [Whole;MCI];

FinalResults = movevars(FinalResults,...
    'Population','Before','Biomarker');

fprintf('\n============================================================\n');
fprintf(' FINAL STEP 47 RESULTS\n');
fprintf('============================================================\n');

disp(FinalResults);

%% ------------------------------------------------------------
% Save
% -------------------------------------------------------------

writetable(FinalResults,...
    fullfile(outputDir,...
    'Step47_LevelVsSlope_Decomposition.csv'));

writetable(WholeSlopeCoefficients,...
    fullfile(outputDir,...
    'Step47_WholeCohort_SlopeCoefficients.csv'));

writetable(MCISlopeCoefficients,...
    fullfile(outputDir,...
    'Step47_MCI_SlopeCoefficients.csv'));

save(fullfile(outputDir,...
    'Step47_LevelVsSlope_Models.mat'),...
    'mdlFull','mdlNoTauSlope','mdlNoPTauSlope',...
    'mdlNoTauMain','mdlNoPTauMain',...
    'mdlFullMCI','mdlNoTauSlopeMCI','mdlNoPTauSlopeMCI',...
    'mdlNoTauMainMCI','mdlNoPTauMainMCI',...
    'FinalResults',...
    'WholeSlopeCoefficients','MCISlopeCoefficients');

%% ------------------------------------------------------------
% Final report
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' STEP 47 COMPLETE\n');
fprintf('============================================================\n');

fprintf('Whole cohort: 143 participants / 379 observations: PASS\n');
fprintf('MCI subgroup: 75 participants / 208 observations: PASS\n');

fprintf('\nPRIMARY SLOPE TESTS\n');

fprintf('\nWHOLE COHORT\n');
fprintf(['Tau PET x Time: LR chi-square(%d) = %.4f, ' ...
         'P = %.6g, Delta AIC = %.4f, Delta BIC = %.4f\n'],...
    df_TauSlope,LR_TauSlope,p_TauSlope,...
    dAIC_TauSlope,dBIC_TauSlope);

fprintf(['p-tau217 x Time: LR chi-square(%d) = %.4f, ' ...
         'P = %.6g, Delta AIC = %.4f, Delta BIC = %.4f\n'],...
    df_PTauSlope,LR_PTauSlope,p_PTauSlope,...
    dAIC_PTauSlope,dBIC_PTauSlope);

fprintf('\nMCI\n');
fprintf(['Tau PET x Time: LR chi-square(%d) = %.4f, ' ...
         'P = %.6g, Delta AIC = %.4f, Delta BIC = %.4f\n'],...
    df_TauSlopeMCI,LR_TauSlopeMCI,p_TauSlopeMCI,...
    dAIC_TauSlopeMCI,dBIC_TauSlopeMCI);

fprintf(['p-tau217 x Time: LR chi-square(%d) = %.4f, ' ...
         'P = %.6g, Delta AIC = %.4f, Delta BIC = %.4f\n'],...
    df_PTauSlopeMCI,LR_PTauSlopeMCI,p_PTauSlopeMCI,...
    dAIC_PTauSlopeMCI,dBIC_PTauSlopeMCI);

fprintf('\nSaved outputs to:\n%s\n',outputDir);

fprintf('\nFiles created:\n');
fprintf('  Step47_LevelVsSlope_Decomposition.csv\n');
fprintf('  Step47_WholeCohort_SlopeCoefficients.csv\n');
fprintf('  Step47_MCI_SlopeCoefficients.csv\n');
fprintf('  Step47_LevelVsSlope_Models.mat\n');

fprintf('============================================================\n');

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function [LR,df,p,dAIC,dBIC] = extractComparison(cmp)

    % Row 1 = reduced model
    % Row 2 = full model

    LR = cmp.LRStat(2);
    df = cmp.deltaDF(2);
    p  = cmp.pValue(2);

    dAIC = cmp.AIC(2) - cmp.AIC(1);
    dBIC = cmp.BIC(2) - cmp.BIC(1);
end

function Out = getSlopeCoefficients(C)

    names = string(C.Name);

    tauIdx = ...
        names=="TimeYears:Tau_z" | ...
        names=="Tau_z:TimeYears";

    ptauIdx = ...
        names=="TimeYears:PTau217_z" | ...
        names=="PTau217_z:TimeYears";

    idx = tauIdx | ptauIdx;

    Out = table(...
        names(idx),...
        C.Estimate(idx),...
        C.SE(idx),...
        C.Lower(idx),...
        C.Upper(idx),...
        C.pValue(idx),...
        'VariableNames',{...
        'Term','Beta','SE','CI_Lower','CI_Upper','PValue'});
end