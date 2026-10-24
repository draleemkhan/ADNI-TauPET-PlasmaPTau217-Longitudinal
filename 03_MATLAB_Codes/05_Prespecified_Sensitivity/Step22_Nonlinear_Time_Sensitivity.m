%% Step22_Nonlinear_Time_Sensitivity.m
% ADNI revised cohort model
% Step 22: Nonlinear-time sensitivity analysis
%
% Purpose:
%   Test whether the primary Tau_z x TimeYears association is robust to
%   allowing nonlinear (quadratic) CDR-SB change over time.
%
% Locked input:
%   D:\output\
%   ADNI_REVISED_COHORT_LOCKED_N143.mat
%
% Models:
%   1) Primary linear-time model (REML)
%   2) Quadratic-time model (REML): adds TimeYears^2
%   3) Tau-curvature model (REML): adds Tau_z x TimeYears^2
%
% Formal fixed-effect comparisons are additionally refitted with ML:
%   Primary vs quadratic-time
%   Quadratic-time vs tau-curvature
%
% Random-effects structure is kept identical throughout:
%   uncorrelated random intercept + random linear time slope by RID.
%
% Outputs are saved to D:\output.

clear; clc;

%% Paths
outDir = 'D:\output';
inFile = fullfile(outDir,'ADNI_REVISED_COHORT_LOCKED_N143.mat');

assert(isfile(inFile), 'Locked cohort file not found: %s', inFile);

S = load(inFile);
assert(isfield(S,'LockedLongitudinal'), ...
    'Expected variable LockedLongitudinal not found in locked MAT file.');

T = S.LockedLongitudinal;

fprintf('\nSTEP 22: NONLINEAR-TIME SENSITIVITY\n');
fprintf('===================================\n');
fprintf('Locked dataset: %d observations, %d participants\n', ...
    height(T), numel(unique(T.RID)));

%% Required variables
req = {'RIDcat','CDRSB','TimeYears','Tau_z','Centiloid_z','PTau217_z', ...
       'Age','SexModel','Education'};
assert(all(ismember(req,T.Properties.VariableNames)), ...
    'One or more required model variables are missing.');

% Verify complete cases for the locked primary model
miss = false(height(T),1);
for k = 1:numel(req)
    v = T.(req{k});
    if isnumeric(v)
        miss = miss | isnan(v);
    elseif iscategorical(v)
        miss = miss | isundefined(v);
    elseif iscell(v)
        miss = miss | cellfun(@isempty,v);
    end
end
assert(~any(miss), 'Unexpected missing values in locked model variables.');

%% Create quadratic time explicitly
T.TimeYears2 = T.TimeYears.^2;

%% Fixed-effect formulas
% Use explicit TimeYears2 variable to keep interpretation transparent.
fixedPrimary = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
                'Age + SexModel + Education + ' ...
                'Tau_z:TimeYears + Centiloid_z:TimeYears + PTau217_z:TimeYears'];

fixedQuad = ['CDRSB ~ TimeYears + TimeYears2 + Tau_z + Centiloid_z + PTau217_z + ' ...
             'Age + SexModel + Education + ' ...
             'Tau_z:TimeYears + Centiloid_z:TimeYears + PTau217_z:TimeYears'];

fixedTauCurve = ['CDRSB ~ TimeYears + TimeYears2 + Tau_z + Centiloid_z + PTau217_z + ' ...
                 'Age + SexModel + Education + ' ...
                 'Tau_z:TimeYears + Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
                 'Tau_z:TimeYears2'];

% Diagonal/uncorrelated random intercept and random linear time slope.
% Separate random terms force zero intercept-slope covariance.
randPart = ' + (1|RIDcat) + (TimeYears-1|RIDcat)';

formPrimary = [fixedPrimary randPart];
formQuad = [fixedQuad randPart];
formTauCurve = [fixedTauCurve randPart];

%% Fit REML sensitivity models
fprintf('\nFitting REML models...\n');

lmePrimary_REML = fitlme(T, formPrimary, ...
    'FitMethod','REML','DummyVarCoding','reference');

lmeQuad_REML = fitlme(T, formQuad, ...
    'FitMethod','REML','DummyVarCoding','reference');

lmeTauCurve_REML = fitlme(T, formTauCurve, ...
    'FitMethod','REML','DummyVarCoding','reference');

%% Helper extraction
getcoef = @(mdl,name) localGetCoef(mdl,name);

p1 = getcoef(lmePrimary_REML,'Tau_z:TimeYears');
p2 = getcoef(lmeQuad_REML,'Tau_z:TimeYears');
q2 = getcoef(lmeQuad_REML,'TimeYears2');
p3 = getcoef(lmeTauCurve_REML,'Tau_z:TimeYears');
q3 = getcoef(lmeTauCurve_REML,'TimeYears2');
tq3 = getcoef(lmeTauCurve_REML,'Tau_z:TimeYears2');

%% REML result table
Sensitivity = table( ...
    ["Primary linear time"; "Quadratic time"; "Quadratic + Tau x Time^2"], ...
    [p1.Estimate; p2.Estimate; p3.Estimate], ...
    [p1.SE; p2.SE; p3.SE], ...
    [p1.Lower; p2.Lower; p3.Lower], ...
    [p1.Upper; p2.Upper; p3.Upper], ...
    [p1.pValue; p2.pValue; p3.pValue], ...
    'VariableNames', {'Model','TauTime_Beta','TauTime_SE','TauTime_CI_Lower', ...
                      'TauTime_CI_Upper','TauTime_P'});

Sensitivity.PercentChangeFromPrimary = ...
    100*(Sensitivity.TauTime_Beta - Sensitivity.TauTime_Beta(1)) ./ ...
    abs(Sensitivity.TauTime_Beta(1));

fprintf('\nTAU x LINEAR-TIME ASSOCIATION (REML)\n');
disp(Sensitivity);

Curvature = table( ...
    ["Quadratic time"; "Quadratic + Tau x Time^2"], ...
    [q2.Estimate; q3.Estimate], ...
    [q2.SE; q3.SE], ...
    [q2.Lower; q3.Lower], ...
    [q2.Upper; q3.Upper], ...
    [q2.pValue; q3.pValue], ...
    [NaN; tq3.Estimate], ...
    [NaN; tq3.SE], ...
    [NaN; tq3.Lower], ...
    [NaN; tq3.Upper], ...
    [NaN; tq3.pValue], ...
    'VariableNames', {'Model','Time2_Beta','Time2_SE','Time2_CI_Lower', ...
    'Time2_CI_Upper','Time2_P','TauTime2_Beta','TauTime2_SE', ...
    'TauTime2_CI_Lower','TauTime2_CI_Upper','TauTime2_P'});

fprintf('\nCURVATURE TERMS (REML)\n');
disp(Curvature);

%% ML refits for valid fixed-effect likelihood-ratio comparisons
fprintf('\nRefitting all three models with ML for nested fixed-effect comparisons...\n');

lmePrimary_ML = fitlme(T, formPrimary, ...
    'FitMethod','ML','DummyVarCoding','reference');

lmeQuad_ML = fitlme(T, formQuad, ...
    'FitMethod','ML','DummyVarCoding','reference');

lmeTauCurve_ML = fitlme(T, formTauCurve, ...
    'FitMethod','ML','DummyVarCoding','reference');

cmpPrimaryVsQuad = compare(lmePrimary_ML,lmeQuad_ML);
cmpQuadVsTauCurve = compare(lmeQuad_ML,lmeTauCurve_ML);

fprintf('\nML LIKELIHOOD-RATIO TEST: PRIMARY vs QUADRATIC TIME\n');
disp(cmpPrimaryVsQuad);

fprintf('\nML LIKELIHOOD-RATIO TEST: QUADRATIC vs QUADRATIC + TAU x TIME^2\n');
disp(cmpQuadVsTauCurve);

%% ML information criteria / log likelihood
models = {lmePrimary_ML,lmeQuad_ML,lmeTauCurve_ML};
names = ["Primary linear time";"Quadratic time";"Quadratic + Tau x Time^2"];

AIC = zeros(3,1); BIC = zeros(3,1); LogLik = zeros(3,1);
for i = 1:3
    AIC(i) = models{i}.ModelCriterion.AIC;
    BIC(i) = models{i}.ModelCriterion.BIC;
    LogLik(i) = models{i}.LogLikelihood;
end

MLFit = table(names,AIC,BIC,LogLik, ...
    'VariableNames',{'Model','AIC','BIC','LogLikelihood'});

fprintf('\nML MODEL FIT\n');
disp(MLFit);

%% Basic convergence / covariance reporting
fprintf('\nRANDOM-EFFECT / RESIDUAL SDs (REML)\n');
fprintf('Primary:\n');
disp(covarianceParameters(lmePrimary_REML));
fprintf('Quadratic:\n');
disp(covarianceParameters(lmeQuad_REML));
fprintf('Quadratic + Tau x Time^2:\n');
disp(covarianceParameters(lmeTauCurve_REML));

%% Save outputs
matFile = fullfile(outDir,'Step22_Nonlinear_Time_Sensitivity_Models.mat');
save(matFile, ...
    'T','Sensitivity','Curvature','MLFit', ...
    'lmePrimary_REML','lmeQuad_REML','lmeTauCurve_REML', ...
    'lmePrimary_ML','lmeQuad_ML','lmeTauCurve_ML', ...
    'cmpPrimaryVsQuad','cmpQuadVsTauCurve', ...
    'formPrimary','formQuad','formTauCurve','-v7.3');

writetable(Sensitivity,fullfile(outDir,'Step22_TauTime_Sensitivity.csv'));
writetable(Curvature,fullfile(outDir,'Step22_Curvature_Terms.csv'));
writetable(MLFit,fullfile(outDir,'Step22_ML_ModelFit.csv'));

fprintf('\nSTEP 22 COMPLETE\n');
fprintf('Saved MAT: %s\n', matFile);
fprintf('Saved CSV summaries to: %s\n', outDir);
fprintf('\nPlease paste the full MATLAB command-window output for interpretation.\n');

%% Local function
function out = localGetCoef(mdl,coefName)
    C = mdl.Coefficients;
    idx = strcmp(C.Name,coefName);

    % MATLAB may reverse interaction-name order.
    if ~any(idx) && contains(coefName,':')
        z = split(coefName,':');
        alt = char(z(2) + ":" + z(1));
        idx = strcmp(C.Name,alt);
    end

    if ~any(idx)
        error('Coefficient "%s" not found. Available names:\n%s', ...
            coefName, strjoin(C.Name,newline));
    end

    r = C(find(idx,1),:);
    out = struct('Estimate',r.Estimate,'SE',r.SE,'Lower',r.Lower, ...
                 'Upper',r.Upper,'pValue',r.pValue);
end
