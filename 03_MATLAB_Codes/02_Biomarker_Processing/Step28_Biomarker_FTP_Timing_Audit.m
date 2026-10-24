%% Step28_Biomarker_FTP_Timing_Audit.m
% ADNI revised cohort
% STEP 28: Audit temporal matching of amyloid PET and Fujirebio p-tau217
% relative to index FTP before deciding on +/-60-day or +/-30-day sensitivity.
%
% This script is intentionally an AUDIT/FEASIBILITY step only.
% It does NOT fit new longitudinal models.
%
% Inputs:
%   D:\output\
%       Step11C_Baseline210_With_Demographics.csv
%   D:\output\
%       ADNI_REVISED_COHORT_LOCKED_N143.mat
%
% Outputs:
%   Step28_Biomarker_FTP_Timing_ParticipantLevel.csv
%   Step28_Biomarker_FTP_Timing_Summary.csv
%   Step28_Biomarker_FTP_Timing_Feasibility.csv
%   Step28_Biomarker_FTP_Timing_Audit.mat
%
% The code searches plausible date/interval variable names rather than silently
% assuming them. If an interval cannot be identified safely, it stops and prints
% all available variable names so the source column can be confirmed.

clear; clc;

outDir = 'D:\output';
baseFile = fullfile(outDir,'Step11C_Baseline210_With_Demographics.csv');
lockedFile = fullfile(outDir,'ADNI_REVISED_COHORT_LOCKED_N143.mat');

assert(isfile(baseFile),'Missing file: %s',baseFile);
assert(isfile(lockedFile),'Missing file: %s',lockedFile);

T210 = readtable(baseFile,'VariableNamingRule','preserve');
S = load(lockedFile,'LockedLongitudinal');
L = S.LockedLongitudinal;

lockedRID = unique(double(L.RID));
assert(numel(lockedRID)==143,'Expected 143 locked participants.');

% Locate RID in baseline table.
ridName = findVar(T210, {'RID','rid'});
assert(~isempty(ridName),'Could not identify RID column in Step11C baseline table.');
rid210 = toNumeric(T210.(ridName));

[tf,loc] = ismember(lockedRID,rid210);
assert(all(tf),'Some locked N=143 RIDs were not found in Step11C baseline table.');
T = T210(loc,:);

fprintf('\n============================================================\n');
fprintf(' STEP 28: BIOMARKER-to-FTP TIMING AUDIT\n');
fprintf('============================================================\n');
fprintf('Locked analytic participants: %d\n',height(T));

%% ---------- FIND INDEX FTP DATE ----------
ftpDateName = findVar(T, { ...
    'IndexFTP_Date','Index FTP Date','FTP_Date','FTP Date', ...
    'IndexFTPDate','TAU_DATE','TauDate','FTPDATE'});

% ---------- DIRECT INTERVAL CANDIDATES ----------
amyIntName = findVar(T, { ...
    'Amyloid_minus_FTP_days','Amyloid_minus_FTP_Days', ...
    'Amyloid_FTP_days','AmyloidFTP_days','AmyloidFTPDays', ...
    'AMYLOID_minus_FTP_days','AmyloidPET_minus_FTP_days', ...
    'AmyloidPET_FTP_Days','AmyloidTimingDays','Amyloid_Days_From_FTP'});

ptauIntName = findVar(T, { ...
    'PTau217_minus_FTP_days','pTau217_minus_FTP_days', ...
    'PTau_minus_FTP_days','PTauFTP_days','PTau217FTPDays', ...
    'Fujirebio_minus_FTP_days','PTauTimingDays', ...
    'PTau217_Days_From_FTP','pT217_minus_FTP_days'});

% ---------- DATE CANDIDATES ----------
amyDateName = findVar(T, { ...
    'Amyloid_Date','Amyloid Date','AmyloidPET_Date','Amyloid PET Date', ...
    'AMYLOID_DATE','AmyloidDate','MatchedAmyloid_Date','Matched Amyloid Date'});

ptauDateName = findVar(T, { ...
    'PTau217_Date','pTau217_Date','PTau_Date','pTau_Date', ...
    'Fujirebio_Date','Fujirebio Date','Plasma_Date','Plasma Date', ...
    'PTau217Date','pT217_Date','pT217Date'});

%% ---------- DERIVE SIGNED INTERVALS ----------
amyDays = nan(height(T),1);
ptauDays = nan(height(T),1);
amySource = "";
ptauSource = "";

if ~isempty(amyIntName)
    amyDays = toNumeric(T.(amyIntName));
    amySource = "direct interval column: " + amyIntName;
elseif ~isempty(ftpDateName) && ~isempty(amyDateName)
    ftpD = toDate(T.(ftpDateName));
    amyD = toDate(T.(amyDateName));
    amyDays = days(amyD-ftpD); % biomarker date minus FTP date
    amySource = "derived as " + amyDateName + " - " + ftpDateName;
end

if ~isempty(ptauIntName)
    ptauDays = toNumeric(T.(ptauIntName));
    ptauSource = "direct interval column: " + ptauIntName;
elseif ~isempty(ftpDateName) && ~isempty(ptauDateName)
    ftpD = toDate(T.(ftpDateName));
    pD = toDate(T.(ptauDateName));
    ptauDays = days(pD-ftpD); % biomarker date minus FTP date
    ptauSource = "derived as " + ptauDateName + " - " + ftpDateName;
end

if all(isnan(amyDays)) || all(isnan(ptauDays))
    fprintf('\nCould not safely identify both biomarker timing intervals.\n');
    fprintf('Available columns in Step11C_Baseline210_With_Demographics.csv:\n\n');
    disp(string(T.Properties.VariableNames)');
    if all(isnan(amyDays))
        fprintf('\nAMYLOID interval not identified.\n');
    end
    if all(isnan(ptauDays))
        fprintf('\nP-TAU217 interval not identified.\n');
    end
    error(['Step 28 stopped rather than guessing. Paste the displayed variable ' ...
           'names so the exact source column(s) can be mapped.']);
end

fprintf('Amyloid timing source: %s\n',amySource);
fprintf('p-tau217 timing source: %s\n',ptauSource);
fprintf('Sign convention: biomarker date minus index FTP date.\n');

%% ---------- VALIDITY CHECK ----------
assert(all(isfinite(amyDays)),'Missing/non-numeric amyloid timing for some N=143 participants.');
assert(all(isfinite(ptauDays)),'Missing/non-numeric p-tau217 timing for some N=143 participants.');

% Cohort was matched within +/-90 days. Stop if interval derivation contradicts that.
if any(abs(amyDays)>90.0001)
    error('Derived amyloid interval exceeds +/-90 days. Check source/sign/date mapping.');
end
if any(abs(ptauDays)>90.0001)
    error('Derived p-tau217 interval exceeds +/-90 days. Check source/sign/date mapping.');
end

%% ---------- SUMMARY ----------
A = summarizeTiming(amyDays,"Amyloid PET");
P = summarizeTiming(ptauDays,"Plasma p-tau217");

Summary = [A;P];

fprintf('\n============================================================\n');
fprintf(' SIGNED TIMING DISTRIBUTIONS (biomarker date - FTP date)\n');
fprintf('============================================================\n');
printTiming(A);
printTiming(P);

%% ---------- WINDOW FEASIBILITY ----------
windows = [90 60 30];
Feasibility = table();

for w = windows
    a = abs(amyDays)<=w;
    p = abs(ptauDays)<=w;
    both = a & p;

    row = table(w, ...
        sum(a),100*mean(a), ...
        sum(p),100*mean(p), ...
        sum(both),100*mean(both), ...
        'VariableNames',{'WindowDays', ...
        'Amyloid_N','Amyloid_Pct', ...
        'PTau217_N','PTau217_Pct', ...
        'Both_N','Both_Pct'});
    Feasibility = [Feasibility;row]; %#ok<AGROW>
end

fprintf('\n============================================================\n');
fprintf(' RESTRICTED-WINDOW FEASIBILITY WITHIN LOCKED N=143\n');
fprintf('============================================================\n');
disp(Feasibility);

%% ---------- DIRECTION COUNTS ----------
Direction = table( ...
    ["Amyloid PET";"Plasma p-tau217"], ...
    [sum(amyDays<0);sum(ptauDays<0)], ...
    [sum(amyDays==0);sum(ptauDays==0)], ...
    [sum(amyDays>0);sum(ptauDays>0)], ...
    'VariableNames',{'Biomarker','BeforeFTP_N','SameDay_N','AfterFTP_N'});

fprintf('\n============================================================\n');
fprintf(' DIRECTION RELATIVE TO INDEX FTP\n');
fprintf('============================================================\n');
disp(Direction);

%% ---------- PARTICIPANT-LEVEL OUTPUT ----------
ParticipantTiming = table(lockedRID,amyDays,ptauDays, ...
    abs(amyDays),abs(ptauDays), ...
    abs(amyDays)<=60,abs(ptauDays)<=60, ...
    abs(amyDays)<=30,abs(ptauDays)<=30, ...
    (abs(amyDays)<=60 & abs(ptauDays)<=60), ...
    (abs(amyDays)<=30 & abs(ptauDays)<=30), ...
    'VariableNames',{'RID','AmyloidMinusFTP_Days','PTau217MinusFTP_Days', ...
    'AbsAmyloidFTP_Days','AbsPTau217FTP_Days', ...
    'AmyloidWithin60','PTau217Within60', ...
    'AmyloidWithin30','PTau217Within30', ...
    'BothWithin60','BothWithin30'});

%% ---------- SAVE ----------
writetable(ParticipantTiming,fullfile(outDir,'Step28_Biomarker_FTP_Timing_ParticipantLevel.csv'));
writetable(Summary,fullfile(outDir,'Step28_Biomarker_FTP_Timing_Summary.csv'));
writetable(Feasibility,fullfile(outDir,'Step28_Biomarker_FTP_Timing_Feasibility.csv'));
writetable(Direction,fullfile(outDir,'Step28_Biomarker_FTP_Timing_Direction.csv'));

save(fullfile(outDir,'Step28_Biomarker_FTP_Timing_Audit.mat'), ...
    'ParticipantTiming','Summary','Feasibility','Direction', ...
    'amyDays','ptauDays','amySource','ptauSource');

fprintf('\n============================================================\n');
fprintf(' STEP 28 AUDIT COMPLETE\n');
fprintf('============================================================\n');
fprintf('Saved outputs to:\n%s\n',outDir);
fprintf('\nPlease paste:\n');
fprintf('1) SIGNED TIMING DISTRIBUTIONS\n');
fprintf('2) RESTRICTED-WINDOW FEASIBILITY\n');
fprintf('3) DIRECTION RELATIVE TO INDEX FTP\n');

%% ================= LOCAL FUNCTIONS =================
function name = findVar(T,candidates)
vars = string(T.Properties.VariableNames);
name = "";
for i=1:numel(candidates)
    j = find(strcmpi(vars,string(candidates{i})),1);
    if ~isempty(j)
        name = vars(j);
        return;
    end
end
% Normalized fallback
normVars = regexprep(lower(vars),'[^a-z0-9]','');
for i=1:numel(candidates)
    c = regexprep(lower(string(candidates{i})),'[^a-z0-9]','');
    j = find(normVars==c,1);
    if ~isempty(j)
        name = vars(j);
        return;
    end
end
end

function x = toNumeric(v)
if isnumeric(v) || islogical(v)
    x = double(v);
elseif iscell(v)
    x = str2double(string(v));
elseif isstring(v) || ischar(v) || iscategorical(v)
    x = str2double(string(v));
else
    x = double(v);
end
x = x(:);
end

function d = toDate(v)
if isdatetime(v)
    d = v;
    return;
end
s = string(v);
formats = {'yyyy-MM-dd','MM/dd/yyyy','dd-MMM-yyyy','yyyyMMdd','dd/MM/yyyy'};
for i=1:numel(formats)
    try
        d = datetime(s,'InputFormat',formats{i});
        if all(~isnat(d))
            return;
        end
    catch
    end
end
try
    d = datetime(s);
catch
    error('Could not parse date column safely.');
end
end

function R = summarizeTiming(x,label)
q = prctile(x,[25 75]);
aq = prctile(abs(x),[25 50 75 90 95]);
R = table(string(label),numel(x), ...
    median(x),q(1),q(2),min(x),max(x), ...
    median(abs(x)),aq(1),aq(3),aq(4),aq(5), ...
    sum(abs(x)<=30),sum(abs(x)<=60),sum(abs(x)<=90), ...
    'VariableNames',{'Biomarker','N','MedianSignedDays','Q1SignedDays', ...
    'Q3SignedDays','MinSignedDays','MaxSignedDays', ...
    'MedianAbsDays','Q1AbsDays','Q3AbsDays','P90AbsDays','P95AbsDays', ...
    'Within30_N','Within60_N','Within90_N'});
end

function printTiming(R)
fprintf('\n%s:\n',R.Biomarker);
fprintf(' N = %d\n',R.N);
fprintf(' Signed interval median %.1f days | IQR %.1f to %.1f | range %.1f to %.1f\n', ...
    R.MedianSignedDays,R.Q1SignedDays,R.Q3SignedDays,R.MinSignedDays,R.MaxSignedDays);
fprintf(' Absolute interval median %.1f days | IQR %.1f to %.1f | P90 %.1f | P95 %.1f\n', ...
    R.MedianAbsDays,R.Q1AbsDays,R.Q3AbsDays,R.P90AbsDays,R.P95AbsDays);
fprintf(' Within +/-30 d: %d/%d (%.1f%%)\n',R.Within30_N,R.N,100*R.Within30_N/R.N);
fprintf(' Within +/-60 d: %d/%d (%.1f%%)\n',R.Within60_N,R.N,100*R.Within60_N/R.N);
fprintf(' Within +/-90 d: %d/%d (%.1f%%)\n',R.Within90_N,R.N,100*R.Within90_N/R.N);
end
