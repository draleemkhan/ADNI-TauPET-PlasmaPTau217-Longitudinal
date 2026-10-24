%% Step27_Baseline_CDRSB_Severity_Adjustment_Sensitivity.m
% ADNI revised cohort
% Step 27: Baseline clinical severity adjustment sensitivity
%
% PURPOSE
% Test whether the primary Tau_z x TimeYears association persists after
% conditioning on participant-level baseline CDR-SB severity.
%
% IMPORTANT ESTIMAND NOTE
% Baseline CDR-SB is already part of the repeated outcome trajectory.
% This sensitivity additionally enters each participant's matched baseline
% CDR-SB as a time-invariant adjustment covariate. It therefore estimates a
% distinct association: tau PET with subsequent CDR-SB trajectory conditional
% on baseline clinical impairment.
%
% PRIMARY SENSITIVITY:
%   Add BaselineCDRSB_z as MAIN EFFECT ONLY.
%
% SECONDARY SLOPE SENSITIVITY:
%   Also add BaselineCDRSB_z x TimeYears, because baseline severity may itself
%   be associated with progression rate. This is reported separately and does
%   not replace the primary sensitivity model.
%
% Random effects:
%   participant-specific random intercept + random slope for LINEAR TimeYears,
%   diagonal covariance, identical to locked primary model.
%
% Estimation:
%   REML for coefficient estimation.
%   ML only for nested model comparisons.

clear; clc;

%% ---------------- PATHS ----------------
outDir = 'D:\output';
lockedFile = fullfile(outDir,'ADNI_REVISED_COHORT_LOCKED_N143.mat');

assert(isfile(lockedFile),'Missing locked longitudinal file: %s',lockedFile);

%% ---------------- LOAD LOCKED DATA ----------------
S = load(lockedFile,'LockedLongitudinal');
L = S.LockedLongitudinal;

assert(height(L)==379,'Expected 379 longitudinal observations.');
assert(numel(unique(double(L.RID)))==143,'Expected 143 participants.');

%% ---------------- IDENTIFY BASELINE OBSERVATION ----------------
% The locked table contains ObservationType. Prefer this rather than assuming
% the observation nearest time zero is baseline.

obsType = lower(strtrim(string(L.ObservationType)));
isBase = contains(obsType,'baseline');

if sum(isBase) ~= 143
    % Conservative fallback: exactly one baseline record per participant,
    % chosen as the record with minimum |TimeDays|. This should not be needed
    % for the locked dataset but prevents silent failure.
    warning(['ObservationType did not identify exactly 143 baseline rows. ' ...
             'Using minimum absolute TimeDays within RID as fallback.']);
    isBase = false(height(L),1);
    rids = unique(double(L.RID));
    for i = 1:numel(rids)
        idx = find(double(L.RID)==rids(i));
        [~,j] = min(abs(double(L.TimeDays(idx))));
        isBase(idx(j)) = true;
    end
end

B = L(isBase,{'RID','CDRSB'});
B.Properties.VariableNames{'CDRSB'} = 'BaselineCDRSB';

assert(height(B)==143,'Expected one baseline CDR-SB value per participant.');
assert(numel(unique(double(B.RID)))==143,'Baseline RID duplication detected.');
assert(all(isfinite(double(B.BaselineCDRSB))),'Missing baseline CDR-SB.');

%% ---------------- STANDARDIZE BASELINE CDR-SB ----------------
% Standardize using the final N=143 participant-level baseline distribution.
% These scaling parameters are fixed for all Step 27 fits.

baselineMean = mean(double(B.BaselineCDRSB));
baselineSD   = std(double(B.BaselineCDRSB),0);

assert(baselineSD>0,'Baseline CDR-SB SD must be >0.');

B.BaselineCDRSB_z = (double(B.BaselineCDRSB)-baselineMean)/baselineSD;

[tf,loc] = ismember(double(L.RID),double(B.RID));
assert(all(tf),'Some longitudinal participants lack baseline CDR-SB.');

L.BaselineCDRSB = double(B.BaselineCDRSB(loc));
L.BaselineCDRSB_z = double(B.BaselineCDRSB_z(loc));

fprintf('\n============================================================\n');
fprintf(' STEP 27: BASELINE CDR-SB SEVERITY ADJUSTMENT\n');
fprintf('============================================================\n');
fprintf('Participants: %d | observations: %d\n',numel(unique(double(L.RID))),height(L));
fprintf('Baseline CDR-SB mean: %.6f | SD: %.6f\n',baselineMean,baselineSD);
fprintf('Baseline CDR-SB median: %.3f | range: %.3f to %.3f\n', ...
    median(B.BaselineCDRSB),min(B.BaselineCDRSB),max(B.BaselineCDRSB));

%% ---------------- MODEL FORMULAS ----------------
fPrimary = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
    'Age + SexModel + Education + Tau_z:TimeYears + ' ...
    'Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
    '(1 + TimeYears|RIDcat)'];

% Main-effect baseline severity adjustment
fBaseMain = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
    'Age + SexModel + Education + BaselineCDRSB_z + ' ...
    'Tau_z:TimeYears + Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
    '(1 + TimeYears|RIDcat)'];

% Secondary model allowing baseline severity to affect progression slope
fBaseSlope = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
    'Age + SexModel + Education + BaselineCDRSB_z + ' ...
    'Tau_z:TimeYears + Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
    'BaselineCDRSB_z:TimeYears + (1 + TimeYears|RIDcat)'];

%% ---------------- REML FITS ----------------
mPrimary = fitlme(L,fPrimary, ...
    'FitMethod','REML','CovariancePattern','Diagonal');

mBaseMain = fitlme(L,fBaseMain, ...
    'FitMethod','REML','CovariancePattern','Diagonal');

mBaseSlope = fitlme(L,fBaseSlope, ...
    'FitMethod','REML','CovariancePattern','Diagonal');

P0 = getTerm(mPrimary,'Tau_z:TimeYears');
P1 = getTerm(mBaseMain,'Tau_z:TimeYears');
P2 = getTerm(mBaseSlope,'Tau_z:TimeYears');

Bmain = getTerm(mBaseMain,'BaselineCDRSB_z');
Bmain2 = getTerm(mBaseSlope,'BaselineCDRSB_z');
Bslope = getTerm(mBaseSlope,'BaselineCDRSB_z:TimeYears');

pct1 = 100*(P1.Estimate-P0.Estimate)/abs(P0.Estimate);
pct2 = 100*(P2.Estimate-P0.Estimate)/abs(P0.Estimate);

fprintf('\n============================================================\n');
fprintf(' TAU x TIME ACROSS BASELINE-SEVERITY MODELS\n');
fprintf('============================================================\n');
fprintf('Primary:\n');
fprintf(' beta %.6f | SE %.6f | 95%% CI %.6f to %.6f | P %.8g\n', ...
    P0.Estimate,P0.SE,P0.Lower,P0.Upper,P0.P);

fprintf('\n+ Baseline CDR-SB main effect:\n');
fprintf(' beta %.6f | SE %.6f | 95%% CI %.6f to %.6f | P %.8g\n', ...
    P1.Estimate,P1.SE,P1.Lower,P1.Upper,P1.P);
fprintf(' Percent change vs primary: %.3f%%\n',pct1);

fprintf('\n+ Baseline CDR-SB main effect AND Baseline CDR-SB x TimeYears:\n');
fprintf(' beta %.6f | SE %.6f | 95%% CI %.6f to %.6f | P %.8g\n', ...
    P2.Estimate,P2.SE,P2.Lower,P2.Upper,P2.P);
fprintf(' Percent change vs primary: %.3f%%\n',pct2);

fprintf('\n============================================================\n');
fprintf(' BASELINE CDR-SB TERMS\n');
fprintf('============================================================\n');
fprintf('Main-effect-only model, BaselineCDRSB_z:\n');
fprintf(' beta %.6f | SE %.6f | 95%% CI %.6f to %.6f | P %.8g\n', ...
    Bmain.Estimate,Bmain.SE,Bmain.Lower,Bmain.Upper,Bmain.P);

fprintf('\nSlope model, BaselineCDRSB_z main effect:\n');
fprintf(' beta %.6f | SE %.6f | 95%% CI %.6f to %.6f | P %.8g\n', ...
    Bmain2.Estimate,Bmain2.SE,Bmain2.Lower,Bmain2.Upper,Bmain2.P);

fprintf('Slope model, BaselineCDRSB_z x TimeYears:\n');
fprintf(' beta %.6f | SE %.6f | 95%% CI %.6f to %.6f | P %.8g\n', ...
    Bslope.Estimate,Bslope.SE,Bslope.Lower,Bslope.Upper,Bslope.P);

%% ---------------- ML NESTED COMPARISONS ----------------
mPrimaryML = fitlme(L,fPrimary, ...
    'FitMethod','ML','CovariancePattern','Diagonal');
mBaseMainML = fitlme(L,fBaseMain, ...
    'FitMethod','ML','CovariancePattern','Diagonal');
mBaseSlopeML = fitlme(L,fBaseSlope, ...
    'FitMethod','ML','CovariancePattern','Diagonal');

LR_main = 2*(mBaseMainML.LogLikelihood-mPrimaryML.LogLikelihood);
p_main = 1-chi2cdf(LR_main,1);

LR_slope = 2*(mBaseSlopeML.LogLikelihood-mBaseMainML.LogLikelihood);
p_slope = 1-chi2cdf(LR_slope,1);

fprintf('\n============================================================\n');
fprintf(' ML NESTED MODEL COMPARISONS\n');
fprintf('============================================================\n');
fprintf('Primary vs + baseline CDR-SB main effect:\n');
fprintf(' LR chi-square(1) = %.6f | P = %.8g\n',LR_main,p_main);
fprintf(' Delta AIC = %.4f | Delta BIC = %.4f\n', ...
    mBaseMainML.ModelCriterion.AIC-mPrimaryML.ModelCriterion.AIC, ...
    mBaseMainML.ModelCriterion.BIC-mPrimaryML.ModelCriterion.BIC);

fprintf('\nBaseline-main model vs + baseline CDR-SB x TimeYears:\n');
fprintf(' LR chi-square(1) = %.6f | P = %.8g\n',LR_slope,p_slope);
fprintf(' Delta AIC = %.4f | Delta BIC = %.4f\n', ...
    mBaseSlopeML.ModelCriterion.AIC-mBaseMainML.ModelCriterion.AIC, ...
    mBaseSlopeML.ModelCriterion.BIC-mBaseMainML.ModelCriterion.BIC);

%% ---------------- SAVE TABLES ----------------
TauResults = table( ...
    ["Primary";"+ Baseline CDR-SB main";"+ Baseline CDR-SB main + slope"], ...
    [P0.Estimate;P1.Estimate;P2.Estimate], ...
    [P0.SE;P1.SE;P2.SE], ...
    [P0.Lower;P1.Lower;P2.Lower], ...
    [P0.Upper;P1.Upper;P2.Upper], ...
    [P0.P;P1.P;P2.P], ...
    [0;pct1;pct2], ...
    'VariableNames',{'Model','TauTimeBeta','SE','CI_Lower','CI_Upper','P','PctChangeVsPrimary'});

BaselineTerms = table( ...
    ["Baseline CDR-SB main (main-only model)"; ...
     "Baseline CDR-SB main (slope model)"; ...
     "Baseline CDR-SB x TimeYears"], ...
    [Bmain.Estimate;Bmain2.Estimate;Bslope.Estimate], ...
    [Bmain.SE;Bmain2.SE;Bslope.SE], ...
    [Bmain.Lower;Bmain2.Lower;Bslope.Lower], ...
    [Bmain.Upper;Bmain2.Upper;Bslope.Upper], ...
    [Bmain.P;Bmain2.P;Bslope.P], ...
    'VariableNames',{'Term','Estimate','SE','CI_Lower','CI_Upper','P'});

MLComparisons = table( ...
    ["Primary vs + baseline CDR-SB main"; ...
     "Baseline-main vs + baseline CDR-SB x Time"], ...
    [LR_main;LR_slope],[1;1],[p_main;p_slope], ...
    [mBaseMainML.ModelCriterion.AIC-mPrimaryML.ModelCriterion.AIC; ...
     mBaseSlopeML.ModelCriterion.AIC-mBaseMainML.ModelCriterion.AIC], ...
    [mBaseMainML.ModelCriterion.BIC-mPrimaryML.ModelCriterion.BIC; ...
     mBaseSlopeML.ModelCriterion.BIC-mBaseMainML.ModelCriterion.BIC], ...
    'VariableNames',{'Comparison','LR_ChiSquare','df','P','DeltaAIC','DeltaBIC'});

writetable(TauResults,fullfile(outDir,'Step27_BaselineCDRSB_TauTime_Results.csv'));
writetable(BaselineTerms,fullfile(outDir,'Step27_BaselineCDRSB_Terms.csv'));
writetable(MLComparisons,fullfile(outDir,'Step27_BaselineCDRSB_ML_Comparisons.csv'));

save(fullfile(outDir,'Step27_Baseline_CDRSB_Severity_Adjustment_Sensitivity.mat'), ...
    'L','B','baselineMean','baselineSD','mPrimary','mBaseMain','mBaseSlope', ...
    'mPrimaryML','mBaseMainML','mBaseSlopeML','TauResults','BaselineTerms', ...
    'MLComparisons','LR_main','p_main','LR_slope','p_slope');

fprintf('\n============================================================\n');
fprintf(' STEP 27 COMPLETE\n');
fprintf('============================================================\n');
fprintf('Saved outputs to:\n%s\n',outDir);
fprintf('\nPlease paste:\n');
fprintf('1) TAU x TIME ACROSS BASELINE-SEVERITY MODELS\n');
fprintf('2) BASELINE CDR-SB TERMS\n');
fprintf('3) ML NESTED MODEL COMPARISONS\n');

%% ---------------- LOCAL FUNCTION ----------------
function R = getTerm(mdl,termName)
C = mdl.Coefficients;
names = string(C.Name);
idx = find(names==string(termName),1);

if isempty(idx) && contains(termName,':')
    p = split(string(termName),':');
    alt = p(2)+":"+p(1);
    idx = find(names==alt,1);
end

assert(~isempty(idx),'Could not find coefficient: %s',termName);

R = struct('Estimate',C.Estimate(idx), ...
           'SE',C.SE(idx), ...
           'Lower',C.Lower(idx), ...
           'Upper',C.Upper(idx), ...
           'P',C.pValue(idx));
end
