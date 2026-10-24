%% Step06_Longitudinal_CDRSB_Construction_v2.m
% ADNI Tau PET revised cohort reconstruction
%
% Robust date-parsing revision of Step06.
%
% Benchmark:
% Participants                         = 143
% Baseline CDR-SB observations         = 143
% Post-index CDR-SB observations       = 236
% Total repeated CDR-SB observations   = 379
% Post-index observations/participant:
%   1 = 73 participants
%   2 = 52 participants
%   3 = 13 participants
%   4 = 5 participants
%
% RULES:
% 1) Index FTP PET defines Time = 0.
% 2) The matched baseline CDR-SB is retained at its ACTUAL assessment time.
% 3) All valid post-index CDR-SB observations from >=3 months through <=4 years are retained.
% 4) No averaging is performed.

clear; clc; close all;

%% Locked path
outputDir = 'D:\output';

finalFile = fullfile(outputDir,'Step05B_Final_143_Participant_Cohort.csv');
postFile  = fullfile(outputDir,'Step05B_All_Qualifying_PostIndex_CDRSB_Long.csv');

assert(isfile(finalFile),'Missing final 143-participant cohort.');
assert(isfile(postFile),'Missing qualifying post-index CDR-SB long file.');

C = readtable(finalFile,'VariableNamingRule','preserve');
Post = readtable(postFile,'VariableNamingRule','preserve');

fprintf('\n====================================================\n');
fprintf(' STEP 6 v2 - LONGITUDINAL CDR-SB CONSTRUCTION\n');
fprintf('====================================================\n');

%% Parse final-cohort dates robustly
if ismember('IndexFTP_Date',C.Properties.VariableNames)
    indexDate = parseDateFlexible(C.IndexFTP_Date);
elseif ismember('SCANDATE',C.Properties.VariableNames)
    indexDate = parseDateFlexible(C.SCANDATE);
else
    error('Index FTP date not found.');
end

if ~ismember('CDRSB_Date',C.Properties.VariableNames)
    error('Matched baseline CDR-SB date not found.');
end
if ~ismember('Baseline_CDRSB',C.Properties.VariableNames)
    error('Baseline_CDRSB not found.');
end

baselineDate = parseDateFlexible(C.CDRSB_Date);

if any(isnat(indexDate))
    error('One or more index FTP dates could not be parsed.');
end
if any(isnat(baselineDate))
    error('One or more baseline CDR-SB dates could not be parsed.');
end

baselineDays = days(baselineDate-indexDate);
baselineYears = baselineDays/365.25;

%% Baseline long-format rows
Base = table();
Base.RID = double(C.RID);
Base.IndexFTP_Date = indexDate;
Base.CDRSB_Date = baselineDate;
Base.TimeDays = baselineDays;
Base.TimeYears = baselineYears;
Base.CDRSB = double(C.Baseline_CDRSB);
Base.ObservationType = repmat("Baseline",height(C),1);

%% Parse post-index long-file dates robustly
if ~ismember('IndexFTP_Date',Post.Properties.VariableNames)
    error('Post-index long file lacks IndexFTP_Date.');
end
if ~ismember('CDRSB_Date',Post.Properties.VariableNames)
    error('Post-index long file lacks CDRSB_Date.');
end

PostLong = table();
PostLong.RID = double(Post.RID);
PostLong.IndexFTP_Date = parseDateFlexible(Post.IndexFTP_Date);
PostLong.CDRSB_Date = parseDateFlexible(Post.CDRSB_Date);

if any(isnat(PostLong.IndexFTP_Date))
    fprintf('WARNING: %d post rows have unparseable IndexFTP_Date.\n',sum(isnat(PostLong.IndexFTP_Date)));
end
if any(isnat(PostLong.CDRSB_Date))
    fprintf('WARNING: %d post rows have unparseable CDRSB_Date.\n',sum(isnat(PostLong.CDRSB_Date)));
end

PostLong.TimeDays = days(PostLong.CDRSB_Date-PostLong.IndexFTP_Date);
PostLong.TimeYears = PostLong.TimeDays/365.25;
PostLong.CDRSB = double(Post.CDRSB);
PostLong.ObservationType = repmat("Post-index",height(Post),1);

%% Verify post-index window
validPost = ~isnat(PostLong.IndexFTP_Date) & ~isnat(PostLong.CDRSB_Date) & ...
            isfinite(PostLong.CDRSB) & ...
            PostLong.TimeYears >= 0.25 & PostLong.TimeYears <= 4.0;

if any(~validPost)
    fprintf('WARNING: %d post-index rows failed date/window validity and were removed.\n', ...
        sum(~validPost));
end
PostLong = PostLong(validPost,:);

%% Combine
Long = [Base; PostLong];
Long = sortrows(Long,{'RID','TimeDays'},{'ascend','ascend'});

%% Core counts
nParticipants = numel(unique(Long.RID));
nBase = height(Base);
nPost = height(PostLong);
nTotal = height(Long);

fprintf('\nCORE COUNTS\n');
fprintf('Participants: %d [benchmark 143]\n',nParticipants);
fprintf('Baseline observations: %d [benchmark 143]\n',nBase);
fprintf('Post-index observations: %d [benchmark 236]\n',nPost);
fprintf('Total repeated observations: %d [benchmark 379]\n',nTotal);

%% Post-index observation counts
rids = unique(Base.RID);
counts = zeros(numel(rids),1);

for i = 1:numel(rids)
    counts(i) = sum(PostLong.RID==rids(i));
end

CountTable = table(rids,counts, ...
    'VariableNames',{'RID','N_PostIndex_CDRSB'});

freqN = (1:4)';
freqParticipants = arrayfun(@(x)sum(counts==x),freqN);

Frequency = table(freqN,freqParticipants, ...
    'VariableNames',{'PostIndex_CDRSB_Observations','Participants'});

fprintf('\nPOST-INDEX CDR-SB OBSERVATIONS PER PARTICIPANT\n');
disp(Frequency);
fprintf('Benchmarks: 1=73 | 2=52 | 3=13 | 4=5\n');

%% Follow-up duration
lastFollow = nan(numel(rids),1);
for i = 1:numel(rids)
    x = PostLong.TimeYears(PostLong.RID==rids(i));
    if ~isempty(x)
        lastFollow(i)=max(x);
    end
end

Followup = table(rids,lastFollow, ...
    'VariableNames',{'RID','LastFollowup_Years'});

medianFU = median(lastFollow,'omitnan');
maxFU = max(lastFollow,[],'omitnan');

fprintf('\nFOLLOW-UP\n');
fprintf('Median last follow-up: %.3f years [benchmark approximately 2.0]\n',medianFU);
fprintf('Maximum follow-up: %.3f years [benchmark <=4.0]\n',maxFU);

%% Baseline timing audit
fprintf('\nBASELINE CDR-SB TIMING RELATIVE TO INDEX FTP\n');
fprintf('Median: %.1f days\n',median(baselineDays,'omitnan'));
fprintf('IQR: %.1f to %.1f days\n',prctile(baselineDays,25),prctile(baselineDays,75));
fprintf('Range: %.0f to %.0f days\n',min(baselineDays),max(baselineDays));
fprintf('Baseline observations retained at actual assessment time.\n');

%% PASS / REVIEW
corePass = nParticipants==143 && nBase==143 && nPost==236 && nTotal==379;
freqPass = isequal(freqParticipants,[73;52;13;5]);
% Mandatory GitHub reproduction checks
assert(corePass, ...
    'Locked longitudinal counts mismatch: expected 143/143/236/379.');
assert(freqPass, ...
    'Locked post-index visit distribution mismatch: expected 73/52/13/5.');


fprintf('\n----------------------------------------------------\n');
if corePass
    fprintf('CORE LONGITUDINAL COUNTS: PASS\n');
else
    fprintf('CORE LONGITUDINAL COUNTS: REVIEW\n');
end

if freqPass
    fprintf('POST-INDEX VISIT DISTRIBUTION: PASS\n');
else
    fprintf('POST-INDEX VISIT DISTRIBUTION: REVIEW\n');
end
fprintf('----------------------------------------------------\n');

%% Save outputs
writetable(Long, ...
    fullfile(outputDir,'Step06_Final_Longitudinal_CDRSB_379_Observations.csv'));

writetable(CountTable, ...
    fullfile(outputDir,'Step06_PostIndex_CDRSB_Count_PerParticipant.csv'));

writetable(Frequency, ...
    fullfile(outputDir,'Step06_PostIndex_CDRSB_Visit_Distribution.csv'));

writetable(Followup, ...
    fullfile(outputDir,'Step06_Followup_Duration_PerParticipant.csv'));

Summary = table(nParticipants,nBase,nPost,nTotal,medianFU,maxFU, ...
    'VariableNames',{'Participants','Baseline_Observations', ...
    'PostIndex_Observations','Total_Observations', ...
    'Median_LastFollowup_Years','Maximum_Followup_Years'});

writetable(Summary, ...
    fullfile(outputDir,'Step06_Longitudinal_CDRSB_Summary.csv'));

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');

%% Local function: robust date parser
function d = parseDateFlexible(x)
    if isdatetime(x)
        d = x;
        return
    end

    s = strtrim(string(x));
    d = NaT(size(s));

    formats = { ...
        'yyyy-MM-dd', ...
        'yyyy-MM-dd HH:mm:ss', ...
        'MM/dd/yyyy', ...
        'dd-MMM-yyyy', ...
        'dd/MM/yyyy'};

    remaining = ~ismissing(s) & strlength(s)>0;

    for k = 1:numel(formats)
        if ~any(remaining), break; end
        try
            tmp = datetime(s(remaining),'InputFormat',formats{k});
            ok = ~isnat(tmp);
            idx = find(remaining);
            d(idx(ok)) = tmp(ok);
            remaining(idx(ok)) = false;
        catch
        end
    end

    % Final MATLAB automatic parsing attempt
    if any(remaining)
        try
            tmp = datetime(s(remaining));
            idx = find(remaining);
            ok = ~isnat(tmp);
            d(idx(ok)) = tmp(ok);
        catch
        end
    end

    d.Format = 'yyyy-MM-dd';
end
