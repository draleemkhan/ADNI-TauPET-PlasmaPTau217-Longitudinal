%% Step04_PTau217_and_Baseline_CDRSB_v2.m
% ADNI Tau PET revised cohort reconstruction
%
% Starting cohort:
%   Step 3B Berkeley amyloid-positive cohort (expected n = 469)
%
% Step 4A:
%   Nearest valid Fujirebio p-tau217 within +/-90 days of index FTP
%   Expected n = 275
%
% Step 4B:
%   Nearest valid CDR-SB within +/-90 days of index FTP
%   Expected n = 211
%
% IMPORTANT:
% - Matching anchor remains INDEX FTP DATE.
% - Plasma field is explicitly pT217_F.
% - Plasma date is explicitly EXAMDATE.
% - CDR-SB field is explicitly CDRSB.
% - CDR date is VISDATE if present, otherwise EXAMDATE.
% - No advanced dementia exclusion yet.
% - No longitudinal follow-up requirement yet.

clear; clc; close all;

%% Locked paths
sourceDir = 'D:\source files';
outputDir = 'D:\output';

if ~exist(outputDir,'dir'), mkdir(outputDir); end

cohortFile = fullfile(outputDir,'Step03B_Amyloid_Positive_Cohort.csv');
ptauFile   = fullfile(sourceDir,'UPENN_PLASMA_FUJIREBIO_QUANTERIX_11Aug2026.csv');
cdrFile    = fullfile(sourceDir,'CDR_11Aug2026.csv');

assert(isfile(cohortFile),'Missing Step03B cohort file.');
assert(isfile(ptauFile),'Missing Fujirebio/Quanterix plasma file.');
assert(isfile(cdrFile),'Missing CDR file.');

C = readtable(cohortFile,'VariableNamingRule','preserve');
P = readtable(ptauFile,'VariableNamingRule','preserve');
D = readtable(cdrFile,'VariableNamingRule','preserve');

fprintf('\n====================================================\n');
fprintf(' STEP 4 v2 - P-TAU217 + BASELINE CDR-SB\n');
fprintf('====================================================\n');
fprintf('Input amyloid-positive cohort: %d\n',height(C));

%% Resolve index FTP date
if ismember('IndexFTP_Date',C.Properties.VariableNames)
    indexDate = datetime(string(C.IndexFTP_Date),'InputFormat','yyyy-MM-dd');
elseif ismember('SCANDATE',C.Properties.VariableNames)
    indexDate = datetime(string(C.SCANDATE),'InputFormat','yyyy-MM-dd');
else
    error('Could not identify index FTP date in Step03B cohort.');
end

%% STEP 4A: Fujirebio p-tau217
requiredP = {'RID','EXAMDATE','pT217_F'};
missingP = setdiff(requiredP,P.Properties.VariableNames);
if ~isempty(missingP)
    error('Plasma file missing required variables: %s',strjoin(missingP,', '));
end

pRID = double(P.RID);
pDate = datetime(string(P.EXAMDATE),'InputFormat','yyyy-MM-dd','Format','yyyy-MM-dd');

if isnumeric(P.pT217_F)
    pValue = double(P.pT217_F);
else
    pValue = str2double(string(P.pT217_F));
end

validP = ~isnan(pRID) & ~isnat(pDate) & isfinite(pValue);

n = height(C);
pMatched = false(n,1);
pRow = nan(n,1);
pGap = nan(n,1);

for i = 1:n
    rid = double(C.RID(i));
    idx = find(validP & pRID==rid);
    if isempty(idx), continue; end

    gaps = days(pDate(idx)-indexDate(i));
    [a,j] = min(abs(gaps));

    if a <= 90
        pMatched(i)=true;
        pRow(i)=idx(j);
        pGap(i)=gaps(j);
    end
end

CP = C(pMatched,:);
sel = pRow(pMatched);

CP.PTau217_Fujirebio = pValue(sel);
CP.PTau217_Date = pDate(sel);
CP.PTau217_minus_FTP_days = pGap(pMatched);
CP.PTau217_FTP_abs_days = abs(pGap(pMatched));

nP = height(CP);

fprintf('\nSTEP 4A: Fujirebio p-tau217 within +/-90 days\n');
fprintf('Matched participants: %d\n',nP);
fprintf('Expected benchmark: 275\n');

if nP==275
    fprintf('STEP 4A STATUS: PASS\n');
else
    fprintf('STEP 4A STATUS: REVIEW (%+d vs benchmark)\n',nP-275);
end

if nP>0
    fprintf('p-tau217 absolute gap: median %.1f days | IQR %.1f-%.1f | max %.0f\n', ...
        median(CP.PTau217_FTP_abs_days,'omitnan'), ...
        prctile(CP.PTau217_FTP_abs_days,25), ...
        prctile(CP.PTau217_FTP_abs_days,75), ...
        max(CP.PTau217_FTP_abs_days));
end

%% STEP 4B: Baseline CDR-SB
if ~ismember('RID',D.Properties.VariableNames)
    error('CDR file missing RID.');
end
if ~ismember('CDRSB',D.Properties.VariableNames)
    error('CDR file missing CDRSB.');
end

if ismember('VISDATE',D.Properties.VariableNames)
    cdrDateVar = 'VISDATE';
elseif ismember('EXAMDATE',D.Properties.VariableNames)
    cdrDateVar = 'EXAMDATE';
else
    error('CDR file has neither VISDATE nor EXAMDATE.');
end

fprintf('\nUsing CDR date variable: %s\n',cdrDateVar);
fprintf('Using CDR-SB variable: CDRSB\n');

dRID = double(D.RID);
dDate = datetime(string(D.(cdrDateVar)), ...
    'InputFormat','yyyy-MM-dd','Format','yyyy-MM-dd');

if isnumeric(D.CDRSB)
    cdrSB = double(D.CDRSB);
else
    cdrSB = str2double(string(D.CDRSB));
end

validD = ~isnan(dRID) & ~isnat(dDate) & isfinite(cdrSB);

n2 = height(CP);
cdrMatched = false(n2,1);
cdrRow = nan(n2,1);
cdrGap = nan(n2,1);

for i = 1:n2
    rid = double(CP.RID(i));

    if ismember('IndexFTP_Date',CP.Properties.VariableNames)
        anchor = datetime(string(CP.IndexFTP_Date(i)),'InputFormat','yyyy-MM-dd');
    else
        anchor = datetime(string(CP.SCANDATE(i)),'InputFormat','yyyy-MM-dd');
    end

    idx = find(validD & dRID==rid);
    if isempty(idx), continue; end

    gaps = days(dDate(idx)-anchor);
    [a,j] = min(abs(gaps));

    if a <= 90
        cdrMatched(i)=true;
        cdrRow(i)=idx(j);
        cdrGap(i)=gaps(j);
    end
end

Final = CP(cdrMatched,:);
selD = cdrRow(cdrMatched);

Final.Baseline_CDRSB = cdrSB(selD);
Final.CDRSB_Date = dDate(selD);
Final.CDRSB_minus_FTP_days = cdrGap(cdrMatched);
Final.CDRSB_FTP_abs_days = abs(cdrGap(cdrMatched));

nCDR = height(Final);

fprintf('\nSTEP 4B: Valid baseline CDR-SB within +/-90 days\n');
fprintf('Matched participants: %d\n',nCDR);
fprintf('Expected benchmark: 211\n');

if nCDR==211
    fprintf('STEP 4B STATUS: PASS\n');
else
    fprintf('STEP 4B STATUS: REVIEW (%+d vs benchmark)\n',nCDR-211);
end

if nCDR>0
    fprintf('CDR-SB absolute gap: median %.1f days | IQR %.1f-%.1f | max %.0f\n', ...
        median(Final.CDRSB_FTP_abs_days,'omitnan'), ...
        prctile(Final.CDRSB_FTP_abs_days,25), ...
        prctile(Final.CDRSB_FTP_abs_days,75), ...
        max(Final.CDRSB_FTP_abs_days));
end

%% Save outputs
writetable(CP,fullfile(outputDir,'Step04A_PTau217_90Day_Cohort.csv'));
writetable(Final,fullfile(outputDir,'Step04B_Baseline_CDRSB_90Day_Cohort.csv'));

NoP = C(~pMatched,:);
writetable(NoP,fullfile(outputDir,'Step04A_No_PTau217_Within90Days.csv'));

NoCDR = CP(~cdrMatched,:);
writetable(NoCDR,fullfile(outputDir,'Step04B_No_CDRSB_Within90Days.csv'));

Flow = table( ...
    ["Amyloid positive"; ...
     "Fujirebio p-tau217 within +/-90 days"; ...
     "Valid baseline CDR-SB within +/-90 days"], ...
    [height(C);nP;nCDR], ...
    [469;275;211], ...
    'VariableNames',{'Step','N_Reproduced','N_Benchmark'});

writetable(Flow,fullfile(outputDir,'Step04_PTau_CDRSB_Flow_Summary.csv'));

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');
