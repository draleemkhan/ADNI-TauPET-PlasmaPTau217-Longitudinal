%% Step05_AdvancedDementia_and_LongitudinalRequirement.m
% ADNI Tau PET revised cohort reconstruction
%
% Starting cohort:
%   Step 4B baseline CDR-SB matched cohort (expected n = 211)
%
% Step 5A:
%   Exclude advanced dementia: baseline CDR-SB >= 9.5
%   Expected n = 210
%
% Step 5B:
%   Require >=1 valid post-index CDR-SB assessment occurring
%   >=3 months and <=4 years after index FTP
%   Expected n = 143
%
% All qualifying repeated post-index CDR-SB observations are retained.
% No averaging is performed.

clear; clc; close all;

%% Locked paths
sourceDir = 'D:\source files';
outputDir = 'D:\output';

if ~exist(outputDir,'dir'), mkdir(outputDir); end

cohortFile = fullfile(outputDir,'Step04B_Baseline_CDRSB_90Day_Cohort.csv');
cdrFile    = fullfile(sourceDir,'CDR_11Aug2026.csv');

assert(isfile(cohortFile),'Missing Step04B cohort file.');
assert(isfile(cdrFile),'Missing CDR source file.');

C = readtable(cohortFile,'VariableNamingRule','preserve');
D = readtable(cdrFile,'VariableNamingRule','preserve');

fprintf('\n====================================================\n');
fprintf(' STEP 5 - ADVANCED DEMENTIA + LONGITUDINAL REQUIREMENT\n');
fprintf('====================================================\n');
fprintf('Input baseline-matched cohort: %d\n',height(C));

%% Resolve index FTP date
if ismember('IndexFTP_Date',C.Properties.VariableNames)
    indexDate = datetime(string(C.IndexFTP_Date),'InputFormat','yyyy-MM-dd');
elseif ismember('SCANDATE',C.Properties.VariableNames)
    indexDate = datetime(string(C.SCANDATE),'InputFormat','yyyy-MM-dd');
else
    error('Could not identify index FTP date.');
end

%% STEP 5A — Advanced dementia exclusion
if ~ismember('Baseline_CDRSB',C.Properties.VariableNames)
    error('Baseline_CDRSB missing from Step04B cohort.');
end

baselineCDR = double(C.Baseline_CDRSB);

keepSeverity = isfinite(baselineCDR) & baselineCDR < 9.5;
C210 = C(keepSeverity,:);

nAfterSeverity = height(C210);
nExcludedSeverity = sum(~keepSeverity);

fprintf('\nSTEP 5A: Exclude baseline CDR-SB >= 9.5\n');
fprintf('Excluded participants: %d\n',nExcludedSeverity);
fprintf('Remaining participants: %d\n',nAfterSeverity);
fprintf('Expected benchmark: 210\n');

if nAfterSeverity == 210
    fprintf('STEP 5A STATUS: PASS\n');
else
    fprintf('STEP 5A STATUS: REVIEW (%+d vs benchmark)\n',nAfterSeverity-210);
end

ExcludedSeverity = C(~keepSeverity,:);

%% Prepare CDR source for longitudinal follow-up
if ~ismember('RID',D.Properties.VariableNames)
    error('CDR source missing RID.');
end
if ~ismember('CDRSB',D.Properties.VariableNames)
    error('CDR source missing CDRSB.');
end

if ismember('VISDATE',D.Properties.VariableNames)
    dateVar = 'VISDATE';
elseif ismember('EXAMDATE',D.Properties.VariableNames)
    dateVar = 'EXAMDATE';
else
    error('CDR source has neither VISDATE nor EXAMDATE.');
end

dRID = double(D.RID);
dDate = datetime(string(D.(dateVar)),'InputFormat','yyyy-MM-dd','Format','yyyy-MM-dd');

if isnumeric(D.CDRSB)
    dCDR = double(D.CDRSB);
else
    dCDR = str2double(string(D.CDRSB));
end

validCDR = ~isnan(dRID) & ~isnat(dDate) & isfinite(dCDR);

%% STEP 5B — Longitudinal eligibility and retention of all qualifying CDR-SB
n = height(C210);

hasFollow = false(n,1);
nQualifying = zeros(n,1);
lastFollowYears = nan(n,1);

% We will build one long row for every qualifying CDR-SB observation
Long = table();

for i = 1:n
    rid = double(C210.RID(i));

    if ismember('IndexFTP_Date',C210.Properties.VariableNames)
        anchor = datetime(string(C210.IndexFTP_Date(i)),'InputFormat','yyyy-MM-dd');
    else
        anchor = datetime(string(C210.SCANDATE(i)),'InputFormat','yyyy-MM-dd');
    end

    idx = find(validCDR & dRID==rid);
    if isempty(idx), continue; end

    deltaDays = days(dDate(idx)-anchor);
    deltaYears = deltaDays / 365.25;

    % Locked qualifying follow-up window:
    % >=3 months approximated as >=0.25 years AND <=4 years
    q = deltaYears >= 0.25 & deltaYears <= 4.0;

    qidx = idx(q);
    qyrs = deltaYears(q);

    if isempty(qidx), continue; end

    hasFollow(i) = true;
    nQualifying(i) = numel(qidx);
    lastFollowYears(i) = max(qyrs);

    tmp = table();
    tmp.RID = repmat(rid,numel(qidx),1);
    tmp.IndexFTP_Date = repmat(anchor,numel(qidx),1);
    tmp.CDRSB_Date = dDate(qidx);
    tmp.TimeYears = qyrs;
    tmp.CDRSB = dCDR(qidx);

    % Preserve selected baseline variables needed downstream
    tmp.Baseline_CDRSB = repmat(double(C210.Baseline_CDRSB(i)),numel(qidx),1);

    if ismember('META_TEMPORAL_SUVR',C210.Properties.VariableNames)
        tmp.TauMetaSUVR = repmat(double(C210.META_TEMPORAL_SUVR(i)),numel(qidx),1);
    end
    if ismember('Centiloid',C210.Properties.VariableNames)
        tmp.Centiloid = repmat(double(C210.Centiloid(i)),numel(qidx),1);
    end
    if ismember('PTau217_Fujirebio',C210.Properties.VariableNames)
        tmp.PTau217_Fujirebio = repmat(double(C210.PTau217_Fujirebio(i)),numel(qidx),1);
    end

    Long = [Long; tmp]; %#ok<AGROW>
end

Final143 = C210(hasFollow,:);
Final143.N_Qualifying_PostIndex_CDRSB = nQualifying(hasFollow);
Final143.Last_Qualifying_Followup_Years = lastFollowYears(hasFollow);

nFinal = height(Final143);
nObs = height(Long);

fprintf('\nSTEP 5B: >=1 post-index CDR-SB from >=3 months to <=4 years\n');
fprintf('Eligible participants: %d\n',nFinal);
fprintf('Expected benchmark: 143\n');
fprintf('Qualifying post-index CDR-SB observations retained: %d\n',nObs);

if nFinal == 143
    fprintf('STEP 5B STATUS: PASS\n');
else
    fprintf('STEP 5B STATUS: REVIEW (%+d vs benchmark)\n',nFinal-143);
end

%% Visit-count audit
if nFinal > 0
    fprintf('\nQualifying post-index CDR-SB count per participant:\n');
    fprintf('Median: %.1f | IQR: %.1f-%.1f | Range: %d-%d\n', ...
        median(Final143.N_Qualifying_PostIndex_CDRSB), ...
        prctile(Final143.N_Qualifying_PostIndex_CDRSB,25), ...
        prctile(Final143.N_Qualifying_PostIndex_CDRSB,75), ...
        min(Final143.N_Qualifying_PostIndex_CDRSB), ...
        max(Final143.N_Qualifying_PostIndex_CDRSB));
end

%% Save outputs
writetable(C210, ...
    fullfile(outputDir,'Step05A_After_AdvancedDementia_Exclusion.csv'));

writetable(ExcludedSeverity, ...
    fullfile(outputDir,'Step05A_Excluded_Baseline_CDRSB_GTE_9p5.csv'));

NoFollow = C210(~hasFollow,:);
writetable(NoFollow, ...
    fullfile(outputDir,'Step05B_No_Qualifying_Longitudinal_CDRSB.csv'));

writetable(Final143, ...
    fullfile(outputDir,'Step05B_Final_143_Participant_Cohort.csv'));

writetable(Long, ...
    fullfile(outputDir,'Step05B_All_Qualifying_PostIndex_CDRSB_Long.csv'));

Flow = table( ...
    ["Baseline CDR-SB within +/-90 days"; ...
     "After baseline CDR-SB <9.5"; ...
     ">=1 post-index CDR-SB >=3 months to <=4 years"], ...
    [height(C);nAfterSeverity;nFinal], ...
    [211;210;143], ...
    'VariableNames',{'Step','N_Reproduced','N_Benchmark'});

writetable(Flow,fullfile(outputDir,'Step05_Cohort_Flow_Summary.csv'));

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');
