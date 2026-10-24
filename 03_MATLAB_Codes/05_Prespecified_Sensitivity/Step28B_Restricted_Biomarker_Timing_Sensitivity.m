%% Step28B_Restricted_Biomarker_Timing_Sensitivity.m
% ADNI revised cohort
% STEP 28B: Restricted biomarker-to-FTP timing sensitivity
%
% Refit the exact locked primary longitudinal mixed model after requiring BOTH:
%   1) amyloid PET, and
%   2) Fujirebio plasma p-tau217
% to be within:
%   A) +/-60 days of index FTP
%   B) +/-30 days of index FTP
%
% IMPORTANT:
% - Uses the existing Tau_z, Centiloid_z and PTau217_z values from the locked
%   N=143 analytic cohort. Z-scores are NOT recalculated after restriction.
%   Thus, a 1-SD biomarker increment has the same scale in all models.
% - Same fixed-effect specification and same diagonal random-intercept/random-
%   linear-time-slope structure as the locked primary model.
% - REML is used for coefficient estimation.
% - No AIC/BIC/LRT comparison is emphasized because restricted models use
%   different participants/observations and are not nested on identical data.

clear; clc;

%% ---------------- PATHS ----------------
outDir = 'D:\output';
lockedFile = fullfile(outDir,'ADNI_REVISED_COHORT_LOCKED_N143.mat');
timingFile = fullfile(outDir,'Step28_Biomarker_FTP_Timing_ParticipantLevel.csv');

assert(isfile(lockedFile),'Missing locked cohort file: %s',lockedFile);
assert(isfile(timingFile),'Missing Step 28 timing file: %s',timingFile);

%% ---------------- LOAD ----------------
S = load(lockedFile,'LockedLongitudinal');
L = S.LockedLongitudinal;
T = readtable(timingFile,'VariableNamingRule','preserve');

assert(height(L)==379,'Expected 379 locked longitudinal observations.');
assert(numel(unique(double(L.RID)))==143,'Expected 143 locked participants.');
assert(height(T)==143,'Expected 143 rows in Step 28 timing file.');

%% ---------------- VERIFY TIMING VARIABLES ----------------
required = {'RID','AmyloidMinusFTP_Days','PTau217MinusFTP_Days'};
for i=1:numel(required)
    assert(ismember(required{i},T.Properties.VariableNames), ...
        'Missing timing variable: %s',required{i});
end

ridTiming = double(T.RID);
amyDays = double(T.AmyloidMinusFTP_Days);
ptauDays = double(T.PTau217MinusFTP_Days);

assert(numel(unique(ridTiming))==143,'Timing table RID duplication detected.');
assert(all(isfinite(amyDays)) && all(isfinite(ptauDays)), ...
    'Timing intervals contain missing/non-finite values.');

keep60RID = ridTiming(abs(amyDays)<=60 & abs(ptauDays)<=60);
keep30RID = ridTiming(abs(amyDays)<=30 & abs(ptauDays)<=30);

assert(numel(keep60RID)==137, ...
    'Expected 137 participants jointly within +/-60 days; found %d.',numel(keep60RID));
assert(numel(keep30RID)==126, ...
    'Expected 126 participants jointly within +/-30 days; found %d.',numel(keep30RID));

L60 = L(ismember(double(L.RID),keep60RID),:);
L30 = L(ismember(double(L.RID),keep30RID),:);

%% ---------------- EXACT PRIMARY MODEL ----------------
formula = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
    'Age + SexModel + Education + Tau_z:TimeYears + ' ...
    'Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
    '(1 + TimeYears|RIDcat)'];

%% ---------------- FIT REML MODELS ----------------
m90 = fitlme(L,formula, ...
    'FitMethod','REML','CovariancePattern','Diagonal');

m60 = fitlme(L60,formula, ...
    'FitMethod','REML','CovariancePattern','Diagonal');

m30 = fitlme(L30,formula, ...
    'FitMethod','REML','CovariancePattern','Diagonal');

R90 = getTerm(m90,'Tau_z:TimeYears');
R60 = getTerm(m60,'Tau_z:TimeYears');
R30 = getTerm(m30,'Tau_z:TimeYears');

pct60 = 100*(R60.Estimate-R90.Estimate)/abs(R90.Estimate);
pct30 = 100*(R30.Estimate-R90.Estimate)/abs(R90.Estimate);

%% ---------------- COUNTS / FOLLOW-UP ----------------
n90 = numel(unique(double(L.RID)));
n60 = numel(unique(double(L60.RID)));
n30 = numel(unique(double(L30.RID)));

o90 = height(L);
o60 = height(L60);
o30 = height(L30);

post90 = sum(~contains(lower(string(L.ObservationType)),'baseline'));
post60 = sum(~contains(lower(string(L60.ObservationType)),'baseline'));
post30 = sum(~contains(lower(string(L30.ObservationType)),'baseline'));

fprintf('\n============================================================\n');
fprintf(' STEP 28B: RESTRICTED BIOMARKER TIMING SENSITIVITY\n');
fprintf('============================================================\n');
fprintf('Scaling: original locked N=143 biomarker z-scores retained.\n');
fprintf('Restriction requires BOTH amyloid PET and p-tau217 within window.\n\n');

fprintf('Primary +/-90 d: N=%d | observations=%d | post-index=%d\n',n90,o90,post90);
fprintf('+/-60 d:        N=%d | observations=%d | post-index=%d\n',n60,o60,post60);
fprintf('+/-30 d:        N=%d | observations=%d | post-index=%d\n',n30,o30,post30);

fprintf('\n============================================================\n');
fprintf(' TAU x TIME RESULTS\n');
fprintf('============================================================\n');

fprintf('Primary +/-90 d:\n');
fprintf(' beta %.6f | SE %.6f | 95%% CI %.6f to %.6f | P %.8g\n', ...
    R90.Estimate,R90.SE,R90.Lower,R90.Upper,R90.P);

fprintf('\nBoth biomarkers within +/-60 d:\n');
fprintf(' beta %.6f | SE %.6f | 95%% CI %.6f to %.6f | P %.8g\n', ...
    R60.Estimate,R60.SE,R60.Lower,R60.Upper,R60.P);
fprintf(' Percent change vs primary: %.3f%%\n',pct60);

fprintf('\nBoth biomarkers within +/-30 d:\n');
fprintf(' beta %.6f | SE %.6f | 95%% CI %.6f to %.6f | P %.8g\n', ...
    R30.Estimate,R30.SE,R30.Lower,R30.Upper,R30.P);
fprintf(' Percent change vs primary: %.3f%%\n',pct30);

%% ---------------- OTHER BIOMARKER x TIME TERMS ----------------
C90 = getTerm(m90,'Centiloid_z:TimeYears');
C60 = getTerm(m60,'Centiloid_z:TimeYears');
C30 = getTerm(m30,'Centiloid_z:TimeYears');

P90 = getTerm(m90,'PTau217_z:TimeYears');
P60 = getTerm(m60,'PTau217_z:TimeYears');
P30 = getTerm(m30,'PTau217_z:TimeYears');

fprintf('\n============================================================\n');
fprintf(' OTHER BIOMARKER x TIME TERMS\n');
fprintf('============================================================\n');
fprintf('Centiloid x Time:\n');
fprintf(' +/-90: beta %.6f | P %.8g\n',C90.Estimate,C90.P);
fprintf(' +/-60: beta %.6f | P %.8g\n',C60.Estimate,C60.P);
fprintf(' +/-30: beta %.6f | P %.8g\n',C30.Estimate,C30.P);

fprintf('\np-tau217 x Time:\n');
fprintf(' +/-90: beta %.6f | P %.8g\n',P90.Estimate,P90.P);
fprintf(' +/-60: beta %.6f | P %.8g\n',P60.Estimate,P60.P);
fprintf(' +/-30: beta %.6f | P %.8g\n',P30.Estimate,P30.P);

%% ---------------- RESULTS TABLE ----------------
Results = table( ...
    ["Primary +/-90 days";"Both biomarkers +/-60 days";"Both biomarkers +/-30 days"], ...
    [n90;n60;n30],[o90;o60;o30],[post90;post60;post30], ...
    [R90.Estimate;R60.Estimate;R30.Estimate], ...
    [R90.SE;R60.SE;R30.SE], ...
    [R90.Lower;R60.Lower;R30.Lower], ...
    [R90.Upper;R60.Upper;R30.Upper], ...
    [R90.P;R60.P;R30.P], ...
    [0;pct60;pct30], ...
    [C90.Estimate;C60.Estimate;C30.Estimate], ...
    [C90.P;C60.P;C30.P], ...
    [P90.Estimate;P60.Estimate;P30.Estimate], ...
    [P90.P;P60.P;P30.P], ...
    'VariableNames',{'Model','Participants','Observations','PostIndexObservations', ...
    'TauTimeBeta','TauTimeSE','TauTimeCI_Lower','TauTimeCI_Upper','TauTimeP', ...
    'TauTimePctChangeVsPrimary','CentiloidTimeBeta','CentiloidTimeP', ...
    'PTau217TimeBeta','PTau217TimeP'});

disp(Results);

%% ---------------- SAVE ----------------
writetable(Results,fullfile(outDir,'Step28B_Restricted_Biomarker_Timing_Sensitivity.csv'));

save(fullfile(outDir,'Step28B_Restricted_Biomarker_Timing_Sensitivity.mat'), ...
    'L60','L30','keep60RID','keep30RID','m90','m60','m30','Results', ...
    'R90','R60','R30','C90','C60','C30','P90','P60','P30');

fprintf('\n============================================================\n');
fprintf(' STEP 28B COMPLETE\n');
fprintf('============================================================\n');
fprintf('Saved outputs to:\n%s\n',outDir);
fprintf('\nPlease paste:\n');
fprintf('1) cohort counts\n');
fprintf('2) TAU x TIME RESULTS\n');
fprintf('3) OTHER BIOMARKER x TIME TERMS\n');

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
