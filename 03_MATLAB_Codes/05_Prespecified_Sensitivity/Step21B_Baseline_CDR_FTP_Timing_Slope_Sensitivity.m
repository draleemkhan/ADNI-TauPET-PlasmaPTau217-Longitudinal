%% Step21B_Baseline_CDR_FTP_Timing_Slope_Sensitivity.m
% ADNI revised cohort
% Step 21B: Baseline CDR-SB-to-FTP timing sensitivity for longitudinal slope
%
% PURPOSE
% Test whether the signed timing interval between baseline CDR-SB and index
% FTP could bias the estimated longitudinal progression slope.
%
% Adds BOTH:
%   1) signed baseline CDR-to-FTP interval as a main effect, and
%   2) signed interval x TimeYears interaction
%
% to the exact locked primary model.
%
% Random effects remain:
%   random intercept + random slope for LINEAR TimeYears only,
%   with diagonal covariance.
%
% Signed interval convention:
%   CDRSB_minus_FTP_days = CDR-SB date minus index FTP date
%   negative = CDR-SB before FTP
%   zero     = same day
%   positive = CDR-SB after FTP
%
% Timing is scaled per 30 days for interpretable coefficients.

clear; clc;

%% Paths
outDir = 'D:\output';
lockedFile = fullfile(outDir,'ADNI_REVISED_COHORT_LOCKED_N143.mat');
baselineFile = fullfile(outDir,'Step11C_Baseline210_With_Demographics.csv');

assert(isfile(lockedFile),'Missing locked longitudinal file.');
assert(isfile(baselineFile),'Missing Step11C baseline file.');

%% Load
S = load(lockedFile,'LockedLongitudinal');
L = S.LockedLongitudinal;

B = readtable(baselineFile,'VariableNamingRule','preserve');

assert(height(L)==379,'Expected 379 longitudinal observations.');
assert(numel(unique(L.RID))==143,'Expected 143 longitudinal participants.');
assert(any(strcmp(B.Properties.VariableNames,'CDRSB_minus_FTP_days')), ...
    'CDRSB_minus_FTP_days not found in Step11C file.');

%% Build participant-level timing table
BT = B(:,{'RID','CDRSB_minus_FTP_days'});
BT = BT(ismember(double(BT.RID),unique(double(L.RID))),:);

assert(height(BT)==143,'Expected timing information for 143 retained participants.');
assert(numel(unique(BT.RID))==143,'RID duplication in baseline timing table.');

BT.TimingDays = double(BT.CDRSB_minus_FTP_days);
BT.Timing30d = BT.TimingDays/30;

assert(all(isfinite(BT.Timing30d)),'Missing/nonfinite baseline timing values.');

%% Map timing to every longitudinal observation
[tf,loc] = ismember(double(L.RID),double(BT.RID));
assert(all(tf),'Some longitudinal participants lack timing information.');

L.TimingDays = BT.TimingDays(loc);
L.Timing30d = BT.Timing30d(loc);

%% Timing QC
participantTiming = BT.TimingDays;

fprintf('\n============================================================\n');
fprintf(' STEP 21B: BASELINE TIMING x TIME SENSITIVITY\n');
fprintf('============================================================\n');
fprintf('Participants: %d\n',height(BT));
fprintf('CDR-SB before FTP: %d (%.1f%%)\n',sum(participantTiming<0),100*mean(participantTiming<0));
fprintf('Same day: %d (%.1f%%)\n',sum(participantTiming==0),100*mean(participantTiming==0));
fprintf('CDR-SB after FTP: %d (%.1f%%)\n',sum(participantTiming>0),100*mean(participantTiming>0));
fprintf('Median signed interval: %.1f days\n',median(participantTiming));
q = prctile(participantTiming,[25 75]);
fprintf('IQR: %.1f to %.1f days\n',q(1),q(2));
fprintf('Range: %.1f to %.1f days\n',min(participantTiming),max(participantTiming));

%% Exact locked primary model
fPrimary = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
    'Age + SexModel + Education + Tau_z:TimeYears + ' ...
    'Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
    '(1 + TimeYears|RIDcat)'];

%% Timing slope-sensitivity model
% Hierarchical with Timing30d main effect and Timing30d x TimeYears.
fTimingSlope = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
    'Age + SexModel + Education + Timing30d + ' ...
    'Tau_z:TimeYears + Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
    'Timing30d:TimeYears + (1 + TimeYears|RIDcat)'];

%% Fit REML diagonal models
mPrimary = fitlme(L,fPrimary, ...
    'FitMethod','REML','CovariancePattern','Diagonal');

mTimingSlope = fitlme(L,fTimingSlope, ...
    'FitMethod','REML','CovariancePattern','Diagonal');

%% Extract terms
P0 = getTerm(mPrimary,'Tau_z:TimeYears');
P1 = getTerm(mTimingSlope,'Tau_z:TimeYears');
Tmain = getTerm(mTimingSlope,'Timing30d');
Tint  = getTerm(mTimingSlope,'Timing30d:TimeYears');

pctChange = 100*(P1.Estimate-P0.Estimate)/abs(P0.Estimate);

fprintf('\n============================================================\n');
fprintf(' PRIMARY vs TIMING-SLOPE-ADJUSTED TAU x TIME\n');
fprintf('============================================================\n');
fprintf('Primary Tau x Time:\n');
fprintf('  beta %.6f | SE %.6f | 95%% CI %.6f to %.6f | P %.8g\n', ...
    P0.Estimate,P0.SE,P0.Lower,P0.Upper,P0.P);

fprintf('Adjusted for Timing30d + Timing30d x TimeYears:\n');
fprintf('  beta %.6f | SE %.6f | 95%% CI %.6f to %.6f | P %.8g\n', ...
    P1.Estimate,P1.SE,P1.Lower,P1.Upper,P1.P);

fprintf('Percent change in Tau x Time beta: %.3f%%\n',pctChange);

fprintf('\nTiming main effect (per 30 days):\n');
fprintf('  beta %.6f | SE %.6f | 95%% CI %.6f to %.6f | P %.8g\n', ...
    Tmain.Estimate,Tmain.SE,Tmain.Lower,Tmain.Upper,Tmain.P);

fprintf('Timing x TimeYears (per 30-day later baseline assessment):\n');
fprintf('  beta %.6f | SE %.6f | 95%% CI %.6f to %.6f | P %.8g\n', ...
    Tint.Estimate,Tint.SE,Tint.Lower,Tint.Upper,Tint.P);

%% Optional ML nested comparison: does adding timing main + slope terms improve fit?
mPrimaryML = fitlme(L,fPrimary, ...
    'FitMethod','ML','CovariancePattern','Diagonal');
mTimingSlopeML = fitlme(L,fTimingSlope, ...
    'FitMethod','ML','CovariancePattern','Diagonal');

LR = 2*(mTimingSlopeML.LogLikelihood-mPrimaryML.LogLikelihood);
df = 2;
pLR = 1-chi2cdf(LR,df);

fprintf('\n============================================================\n');
fprintf(' ML NESTED COMPARISON FOR TWO TIMING TERMS\n');
fprintf('============================================================\n');
fprintf('LR chi-square(%d) = %.6f | P = %.8g\n',df,LR,pLR);
fprintf('Primary ML AIC: %.4f | Timing-adjusted ML AIC: %.4f | Delta AIC: %.4f\n', ...
    mPrimaryML.ModelCriterion.AIC,mTimingSlopeML.ModelCriterion.AIC, ...
    mTimingSlopeML.ModelCriterion.AIC-mPrimaryML.ModelCriterion.AIC);
fprintf('Primary ML BIC: %.4f | Timing-adjusted ML BIC: %.4f | Delta BIC: %.4f\n', ...
    mPrimaryML.ModelCriterion.BIC,mTimingSlopeML.ModelCriterion.BIC, ...
    mTimingSlopeML.ModelCriterion.BIC-mPrimaryML.ModelCriterion.BIC);

%% Save concise result table
Result = table( ...
    ["Primary";"Timing main + timing x time"], ...
    [P0.Estimate;P1.Estimate], ...
    [P0.SE;P1.SE], ...
    [P0.Lower;P1.Lower], ...
    [P0.Upper;P1.Upper], ...
    [P0.P;P1.P], ...
    'VariableNames',{'Model','TauTimeBeta','SE','CI_Lower','CI_Upper','P'});

TimingTerms = table( ...
    ["Timing30d main";"Timing30d x TimeYears"], ...
    [Tmain.Estimate;Tint.Estimate], ...
    [Tmain.SE;Tint.SE], ...
    [Tmain.Lower;Tint.Lower], ...
    [Tmain.Upper;Tint.Upper], ...
    [Tmain.P;Tint.P], ...
    'VariableNames',{'Term','Estimate','SE','CI_Lower','CI_Upper','P'});

writetable(Result,fullfile(outDir,'Step21B_TimingSlope_TauTime_Results.csv'));
writetable(TimingTerms,fullfile(outDir,'Step21B_TimingSlope_TimingTerms.csv'));

save(fullfile(outDir,'Step21B_Baseline_CDR_FTP_Timing_Slope_Sensitivity.mat'), ...
    'L','BT','mPrimary','mTimingSlope','mPrimaryML','mTimingSlopeML', ...
    'Result','TimingTerms','LR','df','pLR','pctChange');

fprintf('\n============================================================\n');
fprintf(' STEP 21B COMPLETE\n');
fprintf('============================================================\n');
fprintf('Saved outputs to:\n%s\n',outDir);
fprintf('\nPlease paste:\n');
fprintf('1) PRIMARY vs TIMING-SLOPE-ADJUSTED TAU x TIME\n');
fprintf('2) Timing main effect and Timing x TimeYears\n');
fprintf('3) ML NESTED COMPARISON FOR TWO TIMING TERMS\n');

%% Local function
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
