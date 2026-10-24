%% Step26A_IPCW_Retention_Sensitivity.m
% ADNI revised cohort
% Step 26A: Inverse-probability-of-retention weighting (IPCW/IPRW) sensitivity
%
% PURPOSE
%   1) Model probability of contributing qualifying longitudinal follow-up
%      among all 210 baseline-eligible participants.
%   2) Derive stabilized inverse-probability-of-retention weights for the
%      143 retained participants.
%   3) Examine positivity / weight stability.
%   4) Refit the locked primary longitudinal model using observation-level
%      participant weights.
%   5) Compare unweighted and weighted Tau_z x TimeYears estimates.
%
% IMPORTANT
%   This addresses selection related to OBSERVED baseline characteristics.
%   It does NOT by itself solve MNAR missingness due to unobserved factors.

clear; clc;

%% ---------------- PATHS ----------------
outDir = 'D:\output';

baselineFile = fullfile(outDir,'Step11C_Baseline210_With_Demographics.csv');
lockedFile   = fullfile(outDir,'ADNI_REVISED_COHORT_LOCKED_N143.mat');

assert(isfile(baselineFile),'Cannot find baseline file: %s',baselineFile);
assert(isfile(lockedFile),'Cannot find locked longitudinal file: %s',lockedFile);

%% ---------------- LOAD DATA ----------------
B = readtable(baselineFile,'VariableNamingRule','preserve');
S = load(lockedFile,'LockedLongitudinal');
L = S.LockedLongitudinal;

fprintf('\n============================================================\n');
fprintf(' STEP 26A: IPCW RETENTION SENSITIVITY\n');
fprintf('============================================================\n');
fprintf('Baseline eligible rows: %d\n',height(B));
fprintf('Locked longitudinal observations: %d\n',height(L));
fprintf('Locked longitudinal participants: %d\n',numel(unique(L.RID)));

assert(height(B)==210,'Expected 210 baseline-eligible participants.');
assert(numel(unique(L.RID))==143,'Expected 143 retained participants.');

%% ---------------- IDENTIFY BASELINE VARIABLES ROBUSTLY ----------------
% This section searches plausible column names because Step11C may use
% slightly different headers. It stops clearly if a required field is absent.

v = B.Properties.VariableNames;

RIDvar = findVar(v, {'RID'});
assert(~isempty(RIDvar),'RID column not found.');

% Retention is defined directly from membership in the locked N=143 cohort.
retainedRID = unique(double(L.RID));
B.Retained = ismember(double(B.(RIDvar)), retainedRID);

fprintf('Retained: %d\n',sum(B.Retained));
fprintf('Non-retained: %d\n',sum(~B.Retained));
assert(sum(B.Retained)==143 && sum(~B.Retained)==67, ...
    'Retention counts do not equal 143/67.');

tauVar = findVar(v, {'META_TEMPORAL_SUVR'});
clVar  = findVar(v, {'Centiloid','CL'});
ptVar  = findVar(v, {'PTau217_Fujirebio','PTau217','pTau217','pT217_F'});
cdrVar = findVar(v, {'CDRSB','Baseline_CDRSB','CDR_SB'});
ageVar = findVar(v, {'Age'});
sexVar = findVar(v, {'Sex','SexModel','Sex_cat'});
eduVar = findVar(v, {'Education','PTEDUCAT','Educ'});
dxVar  = findVar(v, {'Baseline_Diagnosis','Diagnosis','DX'});

required = {tauVar,clVar,ptVar,cdrVar,ageVar,sexVar,eduVar,dxVar};
labels = {'tau','Centiloid','pTau217','baseline CDR-SB','age','sex','education','diagnosis'};
for k=1:numel(required)
    assert(~isempty(required{k}),'Required baseline variable not found: %s',labels{k});
end

fprintf('\nDetected columns:\n');
fprintf('RID: %s\nTau: %s\nCentiloid: %s\npTau217: %s\nCDRSB: %s\nAge: %s\nSex: %s\nEducation: %s\nDiagnosis: %s\n', ...
    RIDvar,tauVar,clVar,ptVar,cdrVar,ageVar,sexVar,eduVar,dxVar);

%% ---------------- FOLLOW-UP OPPORTUNITY ----------------
% Prefer an existing opportunity variable if Step11C contains one.
oppVar = findVar(v, {'FollowupOpportunityYears','FollowUpOpportunityYears', ...
    'OpportunityYears','PotentialFollowupYears','PotentialFollowUpYears', ...
    'Followup_Opportunity_Years','FollowUp_Opportunity_Years'});

if isempty(oppVar)
    % Reconstruct from the previously established administrative cutoff:
    % 07-Aug-2026 minus index FTP date.
    dateVar = findVar(v, {'IndexFTP_Date','Index_FTP_Date','FTP_Date','IndexDate'});
    assert(~isempty(dateVar), ['No follow-up opportunity variable or index FTP date found. ' ...
        'Please paste B.Properties.VariableNames and stop here.']);

    x = B.(dateVar);
    if isdatetime(x)
        indexDate = x;
    elseif isnumeric(x)
        indexDate = datetime(x,'ConvertFrom','datenum');
    else
        indexDate = datetime(string(x));
    end
    adminCutoff = datetime(2026,8,7);
    B.FollowupOpportunityYears = days(adminCutoff-indexDate)/365.25;
    oppVar = 'FollowupOpportunityYears';
    fprintf('Follow-up opportunity reconstructed from index FTP date to 07-Aug-2026.\n');
else
    fprintf('Using existing follow-up opportunity column: %s\n',oppVar);
end

%% ---------------- ANALYSIS DATASET ----------------
A = table;
A.RID = double(B.(RIDvar));
A.Retained = double(B.Retained);
A.Tau = toNumeric(B.(tauVar));
A.Centiloid = toNumeric(B.(clVar));
A.PTau217 = toNumeric(B.(ptVar));
A.CDRSB = toNumeric(B.(cdrVar));
A.Age = toNumeric(B.(ageVar));
A.Education = toNumeric(B.(eduVar));
A.OpportunityYears = toNumeric(B.(oppVar));

A.Sex = categorical(string(B.(sexVar)));
A.Diagnosis = categorical(string(B.(dxVar)));

% Remove unused categorical levels.
A.Sex = removecats(A.Sex);
A.Diagnosis = removecats(A.Diagnosis);

completeMask = all(~ismissing(A(:,{'Tau','Centiloid','PTau217','CDRSB','Age','Education','OpportunityYears'})),2) ...
    & ~ismissing(A.Sex) & ~ismissing(A.Diagnosis);

fprintf('\nRetention-model complete cases: %d/%d\n',sum(completeMask),height(A));
if sum(completeMask) < height(A)
    fprintf('WARNING: %d participants have missing retention-model predictors.\n',sum(~completeMask));
end

AR = A(completeMask,:);

%% ---------------- STANDARDIZE CONTINUOUS PREDICTORS ----------------
% Standardization is across the baseline-eligible analysis set and improves
% numerical conditioning without changing the fitted probabilities.
cont = {'Tau','Centiloid','PTau217','CDRSB','Age','Education','OpportunityYears'};
for k=1:numel(cont)
    x = AR.(cont{k});
    AR.([cont{k} '_z']) = (x-mean(x,'omitnan'))/std(x,0,'omitnan');
end

%% ---------------- RETENTION MODEL ----------------
% Parsimonious prespecified model based only on baseline information and
% administrative follow-up opportunity.
retFormula = ['Retained ~ Tau_z + Centiloid_z + PTau217_z + CDRSB_z + ' ...
              'Age_z + Sex + Education_z + Diagnosis + OpportunityYears_z'];

fprintf('\nRetention model:\n%s\n',retFormula);

glmRet = fitglm(AR,retFormula,'Distribution','binomial','Link','logit');
disp(glmRet.Coefficients);

pRet = predict(glmRet,AR);
pRet = min(max(pRet,1e-6),1-1e-6);

% Stabilized numerator = marginal retention probability in analysis set.
pNum = mean(AR.Retained);

AR.StabilizedWeight = NaN(height(AR),1);
AR.StabilizedWeight(AR.Retained==1) = pNum ./ pRet(AR.Retained==1);

WR = AR(AR.Retained==1,:);

%% ---------------- WEIGHT QC ----------------
w = WR.StabilizedWeight;
q = prctile(w,[0 1 5 25 50 75 95 99 100]);
ESS = (sum(w)^2)/sum(w.^2);

fprintf('\n============================================================\n');
fprintf(' WEIGHT / POSITIVITY QC\n');
fprintf('============================================================\n');
fprintf('Retained participants weighted: %d\n',height(WR));
fprintf('Marginal retention probability: %.6f\n',pNum);
fprintf('Predicted retention probability among retained:\n');
fprintf('  min %.6f | P1 %.6f | P5 %.6f | median %.6f | P95 %.6f | max %.6f\n', ...
    min(pRet(AR.Retained==1)), ...
    prctile(pRet(AR.Retained==1),1), ...
    prctile(pRet(AR.Retained==1),5), ...
    median(pRet(AR.Retained==1)), ...
    prctile(pRet(AR.Retained==1),95), ...
    max(pRet(AR.Retained==1)));

fprintf('Stabilized weights:\n');
fprintf('  min %.4f | P1 %.4f | P5 %.4f | P25 %.4f | median %.4f | P75 %.4f | P95 %.4f | P99 %.4f | max %.4f\n',q);
fprintf('Effective sample size from weights: %.2f of %d\n',ESS,height(WR));

%% ---------------- MAP PARTICIPANT WEIGHTS TO LONGITUDINAL DATA ----------------
[tf,loc] = ismember(double(L.RID), WR.RID);
assert(all(tf),'At least one locked longitudinal RID lacks a retention weight.');

L.IPCW = WR.StabilizedWeight(loc);

% Truncation sensitivity at retained-participant 1st/99th percentiles.
lo99 = prctile(w,1);
hi99 = prctile(w,99);
L.IPCW_Trim99 = min(max(L.IPCW,lo99),hi99);

% More conservative 5th/95th sensitivity.
lo95 = prctile(w,5);
hi95 = prctile(w,95);
L.IPCW_Trim95 = min(max(L.IPCW,lo95),hi95);

%% ---------------- PRIMARY FORMULA ----------------
primaryFormula = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
    'Age + SexModel + Education + Tau_z:TimeYears + ' ...
    'Centiloid_z:TimeYears + PTau217_z:TimeYears + (1 + TimeYears|RIDcat)'];

%% ---------------- FIT MODELS ----------------
% Unweighted exact primary model.
mUnweighted = fitlme(L,primaryFormula, ...
    'FitMethod','REML','CovariancePattern','Diagonal');

% Weighted sensitivity models. MATLAB fitlme accepts observation weights.
mIPCW = fitlme(L,primaryFormula, ...
    'FitMethod','REML','CovariancePattern','Diagonal','Weights',L.IPCW);

mTrim99 = fitlme(L,primaryFormula, ...
    'FitMethod','REML','CovariancePattern','Diagonal','Weights',L.IPCW_Trim99);

mTrim95 = fitlme(L,primaryFormula, ...
    'FitMethod','REML','CovariancePattern','Diagonal','Weights',L.IPCW_Trim95);

%% ---------------- EXTRACT TAU x TIME ----------------
R0 = getTerm(mUnweighted,'Tau_z:TimeYears');
R1 = getTerm(mIPCW,'Tau_z:TimeYears');
R2 = getTerm(mTrim99,'Tau_z:TimeYears');
R3 = getTerm(mTrim95,'Tau_z:TimeYears');

Result = table( ...
    ["Unweighted primary";"IPCW untrimmed";"IPCW 1-99% trimmed";"IPCW 5-95% trimmed"], ...
    [R0.Estimate;R1.Estimate;R2.Estimate;R3.Estimate], ...
    [R0.SE;R1.SE;R2.SE;R3.SE], ...
    [R0.Lower;R1.Lower;R2.Lower;R3.Lower], ...
    [R0.Upper;R1.Upper;R2.Upper;R3.Upper], ...
    [R0.P;R1.P;R2.P;R3.P], ...
    'VariableNames',{'Model','TauTimeBeta','SE','CI_Lower','CI_Upper','P'});

Result.PercentChangeFromPrimary = ...
    100*(Result.TauTimeBeta-Result.TauTimeBeta(1))/abs(Result.TauTimeBeta(1));

fprintf('\n============================================================\n');
fprintf(' TAU x TIME RESULTS\n');
fprintf('============================================================\n');
disp(Result);

%% ---------------- SAVE ----------------
writetable(AR,fullfile(outDir,'Step26A_RetentionModel_All210.csv'));
writetable(WR,fullfile(outDir,'Step26A_IPCW_Retained143.csv'));
writetable(Result,fullfile(outDir,'Step26A_IPCW_PrimaryModel_Results.csv'));

save(fullfile(outDir,'Step26A_IPCW_Retention_Sensitivity.mat'), ...
    'A','AR','WR','L','glmRet','mUnweighted','mIPCW','mTrim99','mTrim95', ...
    'Result','pNum','ESS','q','lo99','hi99','lo95','hi95','retFormula','primaryFormula');

fprintf('\n============================================================\n');
fprintf(' STEP 26A COMPLETE\n');
fprintf('============================================================\n');
fprintf('Saved outputs to:\n%s\n',outDir);
fprintf('\nPlease paste these three sections:\n');
fprintf('1) Detected columns / complete cases\n');
fprintf('2) WEIGHT / POSITIVITY QC\n');
fprintf('3) TAU x TIME RESULTS\n');

%% ================= LOCAL FUNCTIONS =================

function name = findVar(varNames,candidates)
name = '';
vn = string(varNames);
for i=1:numel(candidates)
    hit = find(strcmpi(vn,string(candidates{i})),1);
    if ~isempty(hit)
        name = varNames{hit};
        return;
    end
end
% fallback: normalized exact matching
normVN = lower(regexprep(vn,'[^a-zA-Z0-9]',''));
for i=1:numel(candidates)
    c = lower(regexprep(string(candidates{i}),'[^a-zA-Z0-9]',''));
    hit = find(normVN==c,1);
    if ~isempty(hit)
        name = varNames{hit};
        return;
    end
end
end

function x = toNumeric(x0)
if isnumeric(x0)
    x = double(x0);
elseif islogical(x0)
    x = double(x0);
elseif iscategorical(x0)
    x = str2double(string(x0));
elseif isstring(x0) || iscellstr(x0) || ischar(x0)
    x = str2double(string(x0));
else
    try
        x = double(x0);
    catch
        error('Unable to convert variable to numeric.');
    end
end
x = x(:);
end

function R = getTerm(mdl,termName)
C = mdl.Coefficients;
names = string(C.Name);
idx = find(names==termName,1);
if isempty(idx)
    % Interaction order can occasionally be reversed in display.
    if contains(termName,':')
        p = split(string(termName),':');
        alt = p(2)+":"+p(1);
        idx = find(names==alt,1);
    end
end
assert(~isempty(idx),'Could not find coefficient %s',termName);
R = struct('Estimate',C.Estimate(idx), ...
           'SE',C.SE(idx), ...
           'Lower',C.Lower(idx), ...
           'Upper',C.Upper(idx), ...
           'P',C.pValue(idx));
end
