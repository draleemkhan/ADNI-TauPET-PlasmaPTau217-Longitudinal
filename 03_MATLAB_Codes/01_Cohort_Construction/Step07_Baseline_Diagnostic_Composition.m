%% Step07_Baseline_Diagnostic_Composition.m
% ADNI Tau PET revised cohort reconstruction
%
% Starting cohort:
%   Final 143-participant cohort from Step 5
%
% Benchmark baseline diagnostic composition:
%   Cognitively normal = 33
%   MCI                = 75
%   Dementia           = 35
%   Total              = 143
%
% Baseline diagnosis is matched from DXSUM to the INDEX FTP date.
% The nearest valid diagnostic assessment within +/-90 days is selected.
% This step audits the diagnostic coding before reporting the benchmark.

clear; clc; close all;

%% Locked paths
sourceDir = 'D:\source files';
outputDir = 'D:\output';

cohortFile = fullfile(outputDir,'Step05B_Final_143_Participant_Cohort.csv');
dxFile = fullfile(sourceDir,'DXSUM_11Aug2026.csv');

assert(isfile(cohortFile),'Missing final 143-participant cohort.');
assert(isfile(dxFile),'Missing DXSUM source file.');

C = readtable(cohortFile,'VariableNamingRule','preserve');
D = readtable(dxFile,'VariableNamingRule','preserve');

fprintf('\n====================================================\n');
fprintf(' STEP 7 - BASELINE DIAGNOSTIC COMPOSITION\n');
fprintf('====================================================\n');
fprintf('Final cohort participants: %d\n',height(C));

%% Index FTP date
if ismember('IndexFTP_Date',C.Properties.VariableNames)
    indexDate = parseDateFlexible(C.IndexFTP_Date);
elseif ismember('SCANDATE',C.Properties.VariableNames)
    indexDate = parseDateFlexible(C.SCANDATE);
else
    error('Index FTP date not found.');
end

%% Inspect DXSUM variables
fprintf('\nDXSUM variables:\n');
disp(string(D.Properties.VariableNames)');

assert(ismember('RID',D.Properties.VariableNames),'DXSUM missing RID.');

% Prefer EXAMDATE, otherwise VISDATE
if ismember('EXAMDATE',D.Properties.VariableNames)
    dateVar = 'EXAMDATE';
elseif ismember('VISDATE',D.Properties.VariableNames)
    dateVar = 'VISDATE';
else
    error('DXSUM has neither EXAMDATE nor VISDATE.');
end

%% Identify diagnosis field
vars = string(D.Properties.VariableNames);

preferred = ["DIAGNOSIS","DXCHANGE","DX","DXCURREN"];
dxVar = "";

for k = 1:numel(preferred)
    hit = vars(strcmpi(vars,preferred(k)));
    if ~isempty(hit)
        dxVar = hit(1);
        break
    end
end

if dxVar == ""
    cand = vars(contains(lower(vars),'diag') | startsWith(lower(vars),'dx'));
    fprintf('\nCandidate diagnosis fields:\n');
    disp(cand');
    error('Diagnosis field not automatically resolved. Review candidates above.');
end

fprintf('Using diagnosis date field: %s\n',dateVar);
fprintf('Using diagnosis field: %s\n',dxVar);

%% Convert
dRID = double(D.RID);
dDate = parseDateFlexible(D.(dateVar));
rawDx = D.(char(dxVar));

if isnumeric(rawDx)
    dxNum = double(rawDx);
    dxText = string(dxNum);
else
    dxText = strtrim(string(rawDx));
    dxNum = str2double(dxText);
end

validDx = ~isnan(dRID) & ~isnat(dDate) & ...
          ~ismissing(dxText) & strlength(dxText)>0;

%% Match nearest diagnosis within +/-90 days of index FTP
n = height(C);
matched = false(n,1);
row = nan(n,1);
gap = nan(n,1);

for i = 1:n
    rid = double(C.RID(i));
    idx = find(validDx & dRID==rid);
    if isempty(idx), continue; end

    gaps = days(dDate(idx)-indexDate(i));
    [a,j] = min(abs(gaps));

    if a <= 90
        matched(i)=true;
        row(i)=idx(j);
        gap(i)=gaps(j);
    end
end

fprintf('\nDiagnosis matched within +/-90 days: %d / %d\n',sum(matched),n);

%% Audit raw diagnostic values among matched participants
sel = row(matched);
matchedDxText = dxText(sel);
matchedDxNum = dxNum(sel);

fprintf('\nRAW DIAGNOSIS DISTRIBUTION AMONG MATCHED PARTICIPANTS:\n');
u = unique(matchedDxText);
for k = 1:numel(u)
    fprintf('  %s : %d\n',u(k),sum(matchedDxText==u(k)));
end

%% Map standard ADNI diagnosis codes/text to 3 groups
group = repmat("Unclassified",sum(matched),1);

% Text mappings
txt = upper(strtrim(matchedDxText));

group(contains(txt,"CN") | contains(txt,"NORMAL") | ...
      contains(txt,"NL")) = "Cognitively normal";

group(contains(txt,"MCI")) = "MCI";

group(contains(txt,"DEMENT") | contains(txt,"AD")) = "Dementia";

% Numeric fallback for common DXSUM DIAGNOSIS coding:
% 1 = CN, 2 = MCI, 3 = Dementia
numOK = isfinite(matchedDxNum);
group(numOK & matchedDxNum==1) = "Cognitively normal";
group(numOK & matchedDxNum==2) = "MCI";
group(numOK & matchedDxNum==3) = "Dementia";

%% Attach to cohort
Out = C;
Out.Baseline_Diagnosis = repmat("Unmatched",n,1);
Out.Diagnosis_Date = NaT(n,1);
Out.Diagnosis_minus_FTP_days = nan(n,1);
Out.Raw_Diagnosis = repmat("",n,1);

idxMatched = find(matched);
Out.Baseline_Diagnosis(idxMatched) = group;
Out.Diagnosis_Date(idxMatched) = dDate(sel);
Out.Diagnosis_minus_FTP_days(idxMatched) = gap(matched);
Out.Raw_Diagnosis(idxMatched) = matchedDxText;

%% Counts
nCN = sum(Out.Baseline_Diagnosis=="Cognitively normal");
nMCI = sum(Out.Baseline_Diagnosis=="MCI");
nDem = sum(Out.Baseline_Diagnosis=="Dementia");
nUn = sum(~ismember(Out.Baseline_Diagnosis, ...
    ["Cognitively normal","MCI","Dementia"]));

fprintf('\nBASELINE DIAGNOSTIC COMPOSITION\n');
fprintf('Cognitively normal: %d [benchmark 33]\n',nCN);
fprintf('MCI:                %d [benchmark 75]\n',nMCI);
fprintf('Dementia:           %d [benchmark 35]\n',nDem);
fprintf('Unmatched/unclassified: %d\n',nUn);
fprintf('Total:              %d [benchmark 143]\n',height(Out));

pass = nCN==33 && nMCI==75 && nDem==35 && nUn==0;

fprintf('\n----------------------------------------------------\n');
if pass
    fprintf('BASELINE DIAGNOSTIC COMPOSITION: PASS\n');
else
    fprintf('BASELINE DIAGNOSTIC COMPOSITION: REVIEW\n');
end
fprintf('----------------------------------------------------\n');

%% Save
writetable(Out,fullfile(outputDir,'Step07_Final143_With_Baseline_Diagnosis.csv'));

Composition = table( ...
    ["Cognitively normal";"MCI";"Dementia";"Total"], ...
    [nCN;nMCI;nDem;height(Out)], ...
    [33;75;35;143], ...
    'VariableNames',{'Baseline_Diagnosis','N_Reproduced','N_Benchmark'});

writetable(Composition, ...
    fullfile(outputDir,'Step07_Baseline_Diagnostic_Composition.csv'));

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');

%% Local date parser
function d = parseDateFlexible(x)
    if isdatetime(x)
        d=x; d.Format='yyyy-MM-dd'; return
    end
    s=strtrim(string(x));
    d=NaT(size(s));
    fmts={'yyyy-MM-dd','yyyy-MM-dd HH:mm:ss','MM/dd/yyyy','dd-MMM-yyyy','dd/MM/yyyy'};
    rem=~ismissing(s) & strlength(s)>0;
    for k=1:numel(fmts)
        if ~any(rem), break; end
        try
            tmp=datetime(s(rem),'InputFormat',fmts{k});
            ok=~isnat(tmp); idx=find(rem);
            d(idx(ok))=tmp(ok); rem(idx(ok))=false;
        catch
        end
    end
    if any(rem)
        try
            tmp=datetime(s(rem)); idx=find(rem); ok=~isnat(tmp);
            d(idx(ok))=tmp(ok);
        catch
        end
    end
    d.Format='yyyy-MM-dd';
end
