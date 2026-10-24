%% Step03_Amyloid_Match_and_Positivity.m
% ADNI Tau PET revised cohort reconstruction
%
% Starting cohort:
%   Step 2 earliest usable FTP index cohort (expected n = 1237)
%
% Step 3A:
%   Match the nearest Berkeley amyloid PET to index FTP within +/-90 days
%   Expected n = 1100
%
% Step 3B:
%   Retain Berkeley amyloid-positive classification (AMYLOID_STATUS == 1)
%   Expected n = 469
%
% No plasma p-tau217 or CDR-SB filtering is performed here.

clear; clc; close all;

%% Locked paths
sourceDir = 'D:\source files';
outputDir = 'D:\output';

if ~exist(outputDir,'dir'), mkdir(outputDir); end

indexFile = fullfile(outputDir,'Step02_Earliest_Usable_FTP_Index_Cohort.csv');
amyFile   = fullfile(sourceDir,'All_Subjects_UCBERKELEY_AMY_6MM_12Aug2026.csv');

if ~isfile(indexFile)
    error('Missing Step 2 cohort file: %s',indexFile);
end
if ~isfile(amyFile)
    error('Missing Berkeley amyloid file: %s',amyFile);
end

Index = readtable(indexFile,'VariableNamingRule','preserve');
Amy   = readtable(amyFile,'VariableNamingRule','preserve');

%% Verify required fields
reqI = {'RID','SCANDATE','META_TEMPORAL_SUVR'};
reqA = {'RID','SCANDATE','TRACER','AMYLOID_STATUS','CENTILOIDS'};

mi = setdiff(reqI,Index.Properties.VariableNames);
ma = setdiff(reqA,Amy.Properties.VariableNames);

if ~isempty(mi)
    error('Index cohort missing variables: %s',strjoin(mi,', '));
end
if ~isempty(ma)
    error('Amyloid file missing variables: %s',strjoin(ma,', '));
end

%% Parse dates
Index.IndexFTP_Date = datetime(string(Index.SCANDATE), ...
    'InputFormat','yyyy-MM-dd','Format','yyyy-MM-dd');

Amy.Amy_Date = datetime(string(Amy.SCANDATE), ...
    'InputFormat','yyyy-MM-dd','Format','yyyy-MM-dd');

% Basic valid amyloid records for temporal matching
amyRID = double(Amy.RID);
validAmy = ~isnan(amyRID) & ~isnat(Amy.Amy_Date);

AmyValid = Amy(validAmy,:);

%% Match closest amyloid PET within +/-90 days for each index FTP RID
nIndex = height(Index);

matched = false(nIndex,1);
amyRow = nan(nIndex,1);
signedGap = nan(nIndex,1);
absGap = nan(nIndex,1);

for i = 1:nIndex
    rid = double(Index.RID(i));

    idx = find(double(AmyValid.RID)==rid);

    if isempty(idx)
        continue
    end

    gaps = days(AmyValid.Amy_Date(idx) - Index.IndexFTP_Date(i));
    [minAbs,j] = min(abs(gaps));

    if minAbs <= 90
        matched(i) = true;
        amyRow(i) = idx(j);
        signedGap(i) = gaps(j);   % + = amyloid after FTP
        absGap(i) = minAbs;
    end
end

%% Build matched cohort
M = Index(matched,:);

sel = amyRow(matched);
A = AmyValid(sel,:);

% Append only explicitly named amyloid variables
M.Amyloid_Date = A.Amy_Date;
M.Amyloid_Tracer = string(A.TRACER);
M.Amyloid_Status_Berkeley = double(A.AMYLOID_STATUS);
M.Centiloid = double(A.CENTILOIDS);
M.Amyloid_minus_FTP_days = signedGap(matched);
M.Amyloid_FTP_abs_days = absGap(matched);

nMatched = height(M);

fprintf('\n====================================================\n');
fprintf(' STEP 3 - AMYLOID MATCHING AND POSITIVITY\n');
fprintf('====================================================\n');
fprintf('Input earliest usable FTP cohort: %d\n',nIndex);

fprintf('\nSTEP 3A: Amyloid PET within +/-90 days\n');
fprintf('Matched participants: %d\n',nMatched);
fprintf('Expected benchmark: 1100\n');

if nMatched == 1100
    fprintf('STEP 3A STATUS: PASS\n');
else
    fprintf('STEP 3A STATUS: REVIEW (%+d vs benchmark)\n',nMatched-1100);
end

%% Gap audit
if nMatched > 0
    fprintf('Amyloid-FTP absolute gap: median %.1f days | IQR %.1f-%.1f | max %.0f\n', ...
        median(M.Amyloid_FTP_abs_days,'omitnan'), ...
        prctile(M.Amyloid_FTP_abs_days,25), ...
        prctile(M.Amyloid_FTP_abs_days,75), ...
        max(M.Amyloid_FTP_abs_days));
end

%% Berkeley amyloid positivity
% Source file provides AMYLOID_STATUS and AMYLOID_STATUS_COMPOSITE_REF.
% Primary cohort uses Berkeley AMYLOID_STATUS == 1.
validStatus = isfinite(M.Amyloid_Status_Berkeley);
AmyPos = M(validStatus & M.Amyloid_Status_Berkeley==1,:);

nValidStatus = sum(validStatus);
nPositive = height(AmyPos);

fprintf('\nSTEP 3B: Berkeley amyloid-positive classification\n');
fprintf('Matched participants with valid AMYLOID_STATUS: %d\n',nValidStatus);
fprintf('Amyloid-positive participants: %d\n',nPositive);
fprintf('Expected benchmark: 469\n');

if nPositive == 469
    fprintf('STEP 3B STATUS: PASS\n');
else
    fprintf('STEP 3B STATUS: REVIEW (%+d vs benchmark)\n',nPositive-469);
end

%% Amyloid status distribution
statusVals = M.Amyloid_Status_Berkeley;
fprintf('\nBerkeley AMYLOID_STATUS distribution among matched participants:\n');
u = unique(statusVals(isfinite(statusVals)));
for k = 1:numel(u)
    fprintf('  Status %.0f: %d\n',u(k),sum(statusVals==u(k)));
end
fprintf('  Missing/non-numeric: %d\n',sum(~isfinite(statusVals)));

%% Save outputs
writetable(M,fullfile(outputDir,'Step03A_Amyloid_Matched_90Day_Cohort.csv'));
writetable(AmyPos,fullfile(outputDir,'Step03B_Amyloid_Positive_Cohort.csv'));

Unmatched = Index(~matched,:);
writetable(Unmatched,fullfile(outputDir,'Step03A_No_Amyloid_Within_90Days.csv'));

Flow = table( ...
    ["Earliest usable FTP index"; ...
     "Amyloid PET within +/-90 days"; ...
     "Berkeley amyloid positive"], ...
    [nIndex;nMatched;nPositive], ...
    [1237;1100;469], ...
    'VariableNames',{'Step','N_Reproduced','N_Benchmark'});

writetable(Flow,fullfile(outputDir,'Step03_Amyloid_Flow_Summary.csv'));

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');
