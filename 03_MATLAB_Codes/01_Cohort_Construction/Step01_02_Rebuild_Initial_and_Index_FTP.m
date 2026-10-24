%% Step01_02_Rebuild_Initial_and_Index_FTP.m
% ADNI Tau PET revised cohort reconstruction
% Reproduces:
%   Step 1: >=1 flortaucipir (FTP) PET participant
%   Step 2: earliest usable FTP index scan
%
% Authoritative tau source: NON-PVC UCBERKELEY_TAU_6MM
% No amyloid, plasma, CDR-SB, diagnosis, or longitudinal filtering here.

clear; clc; close all;

%% Locked paths
sourceDir = 'D:\source files';
outputDir = 'D:\output';

if ~exist(outputDir,'dir'), mkdir(outputDir); end

%% Locate authoritative non-PVC tau file
candidates = dir(fullfile(sourceDir,'*UCBERKELEY_TAU_6MM*.csv'));
names = string({candidates.name});

% Exclude TAUPVC files explicitly
keep = ~contains(upper(names),'TAUPVC');
candidates = candidates(keep);

if isempty(candidates)
    error('No non-PVC UCBERKELEY_TAU_6MM CSV found in sourceDir.');
elseif numel(candidates) > 1
    fprintf('Multiple non-PVC tau files found:\n');
    disp(string({candidates.name})');
    error('Please retain/identify one authoritative non-PVC tau file.');
end

tauFile = fullfile(sourceDir,candidates(1).name);
fprintf('\nUsing authoritative tau file:\n%s\n',tauFile);

Tau = readtable(tauFile,'VariableNamingRule','preserve');

%% Required fields
req = {'RID','SCANDATE','TRACER','META_TEMPORAL_SUVR'};
missing = setdiff(req,Tau.Properties.VariableNames);
if ~isempty(missing)
    error('Required tau variables missing: %s',strjoin(missing,', '));
end

%% STEP 1 — all participants with >=1 FTP PET
tracer = upper(strtrim(string(Tau.TRACER)));
isFTP = tracer == "FTP";

FTP = Tau(isFTP,:);

rid = double(FTP.RID);
validRID = ~isnan(rid);
initialRIDs = unique(rid(validRID));

nInitial = numel(initialRIDs);

fprintf('\n====================================================\n');
fprintf(' REVISED COHORT: STEP 1 + STEP 2 FTP REPRODUCTION\n');
fprintf('====================================================\n');
fprintf('Total rows in tau source: %d\n',height(Tau));
fprintf('FTP rows: %d\n',height(FTP));
fprintf('\nSTEP 1: Participants with >=1 FTP PET: %d\n',nInitial);
fprintf('Expected benchmark: 1246\n');

if nInitial == 1246
    fprintf('STEP 1 STATUS: PASS\n');
else
    fprintf('STEP 1 STATUS: REVIEW (%+d vs benchmark)\n',nInitial-1246);
end

%% STEP 2 — define usable quantitative FTP
scanDate = datetime(string(FTP.SCANDATE),'InputFormat','yyyy-MM-dd', ...
    'Format','yyyy-MM-dd');

suvr = double(FTP.META_TEMPORAL_SUVR);

usable = validRID & ~isnat(scanDate) & isfinite(suvr);

FTPusable = FTP(usable,:);
FTPusable.ScanDate_dt = scanDate(usable);

usableRIDs = unique(double(FTPusable.RID));
nUsableParticipants = numel(usableRIDs);

% Earliest usable scan per RID
FTPusable = sortrows(FTPusable,{'RID','ScanDate_dt'},{'ascend','ascend'});
[~,ia] = unique(double(FTPusable.RID),'stable');
IndexFTP = FTPusable(ia,:);

nIndex = height(IndexFTP);

fprintf('\nSTEP 2: Participants with >=1 usable FTP result: %d\n', ...
    nUsableParticipants);
fprintf('Earliest usable FTP index participants: %d\n',nIndex);
fprintf('Expected benchmark: 1237\n');

if nIndex == 1237
    fprintf('STEP 2 STATUS: PASS\n');
else
    fprintf('STEP 2 STATUS: REVIEW (%+d vs benchmark)\n',nIndex-1237);
end

fprintf('\nExcluded between Step 1 and Step 2: %d participants\n', ...
    nInitial-nIndex);

%% Participant-level exclusion audit
excludedRIDs = setdiff(initialRIDs,usableRIDs);
Excluded = table(excludedRIDs,'VariableNames',{'RID'});

%% Save initial participant list
Initial = table(initialRIDs,'VariableNames',{'RID'});
writetable(Initial,fullfile(outputDir,'Step01_Initial_FTP_Participants.csv'));

%% Save index cohort
preferred = {'RID','PTID','VISCODE','VISCODE2','SCANDATE','TRACER', ...
             'META_TEMPORAL_SUVR','ScanDate_dt'};
keepVars = intersect(preferred,IndexFTP.Properties.VariableNames,'stable');
IndexCompact = IndexFTP(:,keepVars);

writetable(IndexCompact, ...
    fullfile(outputDir,'Step02_Earliest_Usable_FTP_Index_Cohort.csv'));

writetable(Excluded, ...
    fullfile(outputDir,'Step02_Excluded_NoUsable_FTP.csv'));

%% Save flow summary
Flow = table( ...
    ["Initial FTP participants";"Earliest usable FTP index"], ...
    [nInitial;nIndex], ...
    [1246;1237], ...
    'VariableNames',{'Step','N_Reproduced','N_Benchmark'});

writetable(Flow,fullfile(outputDir,'Step01_02_FTP_Cohort_Flow.csv'));

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');
