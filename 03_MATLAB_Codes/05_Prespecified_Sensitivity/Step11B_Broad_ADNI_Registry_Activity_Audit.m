%% Step11B_Broad_ADNI_Registry_Activity_Audit.m
% ADNI Tau PET revised cohort analysis
%
% PURPOSE:
% Reconstruct "later registry activity" among the 67 baseline-eligible
% participants who lacked qualifying longitudinal CDR-SB, using ALL 14
% available ADNI source CSV files rather than CDR/DXSUM alone.
%
% This is an audit only. It does not alter the locked N=143 cohort.
%
% Historical benchmark:
%   N baseline eligible = 210
%   N longitudinal eligible = 143
%   N without qualifying longitudinal CDR-SB = 67
%   Verified later registry activity = 22
%   Verified no later registry activity = 45
%   Verified CN among later-registry group = 18/22
%
% A source file contributes registry evidence only if:
%   1) it contains RID;
%   2) it contains a parseable assessment/scan/date field;
%   3) the dated record occurs >=3 months and <=4 years after index FTP.
%
% The script saves participant-level and source-level audit tables.

clear; clc; close all;

sourceDir = 'D:\source files';
outputDir = 'D:\output';

baseFile = fullfile(outputDir,'Step05A_After_AdvancedDementia_Exclusion.csv');
finalFile = fullfile(outputDir,'Step08_Final143_With_Demographics.csv');

assert(isfile(baseFile),'Missing Step05A baseline N=210 cohort.');
assert(isfile(finalFile),'Missing Step08 final N=143 cohort.');

B = readtable(baseFile,'VariableNamingRule','preserve');
F = readtable(finalFile,'VariableNamingRule','preserve');

if ismember('IndexFTP_Date',B.Properties.VariableNames)
    t0 = pdate(B.IndexFTP_Date);
elseif ismember('SCANDATE',B.Properties.VariableNames)
    t0 = pdate(B.SCANDATE);
else
    error('Index FTP date field not found.');
end

B.IndexFTP_Date_Audit = t0;
B.LongitudinalEligible = ismember(double(B.RID),double(F.RID));
nonret = ~B.LongitudinalEligible;

fprintf('\n====================================================\n');
fprintf(' STEP 11B - BROAD ADNI REGISTRY ACTIVITY AUDIT\n');
fprintf('====================================================\n');
fprintf('Baseline eligible: %d [210]\n',height(B));
fprintf('Longitudinal eligible: %d [143]\n',sum(B.LongitudinalEligible));
fprintf('Non-retained: %d [67]\n',sum(nonret));

%% Deterministic source selection — locked 14-file inventory
inventoryFile = fullfile(fileparts(mfilename('fullpath')), ...
    '..','..','01_Documentation','Data_Access', ...
    'ADNI_Source_Data_Inventory_14.csv');

assert(isfile(inventoryFile),'Locked source inventory missing.');
I = readtable(inventoryFile,'TextType','string');
assert(ismember('FileName',I.Properties.VariableNames));
sourceNames = string(I.FileName);
assert(numel(sourceNames)==14);
assert(numel(unique(sourceNames))==14);

csvs = struct('folder',{},'name',{});
for k = 1:numel(sourceNames)
    assert(isfile(fullfile(sourceDir,sourceNames(k))), ...
        'Required ADNI source CSV missing: %s',sourceNames(k));
    csvs(k).folder = sourceDir;
    csvs(k).name = char(sourceNames(k));
end

% Participant x source evidence matrix
srcNames = strings(numel(csvs),1);
srcEvidence = false(height(B),numel(csvs));
SourceAudit = table();

for s = 1:numel(csvs)
    fn = fullfile(csvs(s).folder,csvs(s).name);
    srcNames(s) = string(csvs(s).name);

    try
        T = readtable(fn,'VariableNamingRule','preserve');
    catch ME
        row = table(string(csvs(s).name),"READ_ERROR","",0,0, ...
            'VariableNames',{'File','Status','DateField','Rows','NonRetainedWithLaterRecord'});
        SourceAudit = [SourceAudit;row]; %#ok<AGROW>
        fprintf('[READ ERROR] %s : %s\n',csvs(s).name,ME.message);
        continue
    end

    vars = string(T.Properties.VariableNames);

    if ~ismember("RID",vars)
        row = table(string(csvs(s).name),"NO_RID","",height(T),0, ...
            'VariableNames',{'File','Status','DateField','Rows','NonRetainedWithLaterRecord'});
        SourceAudit = [SourceAudit;row]; %#ok<AGROW>
        fprintf('[SKIP] %-55s no RID\n',csvs(s).name);
        continue
    end

    % Candidate date fields in priority order.
    priorities = ["EXAMDATE","SCANDATE","VISDATE","USERDATE","USERDATE2", ...
                  "COLPROT","RUNDATE","DATE","StudyDate","IMAGE_DATE","image_date"];

    dateField = "";
    for p = priorities
        idx = find(strcmpi(vars,p),1);
        if ~isempty(idx)
            candidate = vars(idx);
            dtest = pdate(T.(char(candidate)));
            if sum(~isnat(dtest)) > 0
                dateField = candidate;
                break
            end
        end
    end

    % If priority names fail, inspect any field containing DATE.
    if dateField==""
        dateVars = vars(contains(upper(vars),"DATE"));
        bestN = 0;
        for dv = dateVars
            dtest = pdate(T.(char(dv)));
            n = sum(~isnat(dtest));
            if n > bestN
                bestN = n;
                dateField = dv;
            end
        end
    end

    if dateField==""
        row = table(string(csvs(s).name),"NO_PARSEABLE_DATE","",height(T),0, ...
            'VariableNames',{'File','Status','DateField','Rows','NonRetainedWithLaterRecord'});
        SourceAudit = [SourceAudit;row]; %#ok<AGROW>
        fprintf('[SKIP] %-55s no parseable date\n',csvs(s).name);
        continue
    end

    rid = toNum(T.RID);
    d = pdate(T.(char(dateField)));

    for i = find(nonret)'
        q = rid==double(B.RID(i)) & ~isnat(d);
        if any(q)
            dt = days(d(q)-t0(i))/365.25;
            srcEvidence(i,s) = any(dt>=0.25 & dt<=4);
        end
    end

    nEvidence = sum(nonret & srcEvidence(:,s));

    row = table(string(csvs(s).name),"USED",string(dateField),height(T),nEvidence, ...
        'VariableNames',{'File','Status','DateField','Rows','NonRetainedWithLaterRecord'});
    SourceAudit = [SourceAudit;row]; %#ok<AGROW>

    fprintf('[USED] %-55s date=%-12s non-retained later=%d\n', ...
        csvs(s).name,dateField,nEvidence);
end

%% Broad registry flag
broadLater = any(srcEvidence,2);
B.BroadLaterRegistryActivity = broadLater;

nLater = sum(nonret & broadLater);
nNoLater = sum(nonret & ~broadLater);

fprintf('\n----------------------------------------------------\n');
fprintf('BROAD ALL-SOURCE RESULT\n');
fprintf('Later registry activity, no qualifying CDR-SB: %d [verified 22]\n',nLater);
fprintf('No later registry activity: %d [verified 45]\n',nNoLater);

%% Baseline diagnosis nearest +/-90 d
dxFile=fullfile(sourceDir,'DXSUM_11Aug2026.csv');
DX=readtable(dxFile,'VariableNamingRule','preserve');
drid=toNum(DX.RID);
dd=pdate(DX.EXAMDATE);
dc=toNum(DX.DIAGNOSIS);

bdx=nan(height(B),1);
for i=1:height(B)
    q=find(drid==double(B.RID(i)) & ~isnat(dd) & isfinite(dc));
    if isempty(q), continue; end
    gap=days(dd(q)-t0(i));
    [m,j]=min(abs(gap));
    if m<=90, bdx(i)=dc(q(j)); end
end

B.BaselineDiagnosisCode=bdx;
B.BaselineDiagnosis=repmat("Unclassified",height(B),1);
B.BaselineDiagnosis(bdx==1)="Cognitively normal";
B.BaselineDiagnosis(bdx==2)="MCI";
B.BaselineDiagnosis(bdx==3)="Dementia";

nCN=sum(nonret & broadLater & B.BaselineDiagnosis=="Cognitively normal");
fprintf('CN among later-registry/non-retained: %d / %d [verified 18/22]\n',nCN,nLater);

%% Incremental contribution by source
% Show which files add participants beyond the original CDR + DXSUM definition.
cdrCol=find(strcmpi(srcNames,'CDR_11Aug2026.csv'),1);
dxCol=find(strcmpi(srcNames,'DXSUM_11Aug2026.csv'),1);

core=false(height(B),1);
if ~isempty(cdrCol), core=core|srcEvidence(:,cdrCol); end
if ~isempty(dxCol), core=core|srcEvidence(:,dxCol); end

fprintf('\nCORE CDR/DXSUM later activity: %d\n',sum(nonret&core));
fprintf('Additional participants detected by other sources: %d\n', ...
    sum(nonret & broadLater & ~core));

AdditionalBySource=table();
for s=1:numel(csvs)
    if s==cdrCol || s==dxCol, continue; end
    add = nonret & srcEvidence(:,s) & ~core;
    if any(add)
        row=table(srcNames(s),sum(add), ...
            'VariableNames',{'SourceFile','AdditionalBeyond_CDR_DXSUM'});
        AdditionalBySource=[AdditionalBySource;row]; %#ok<AGROW>
    end
end

if ~isempty(AdditionalBySource)
    fprintf('\nSOURCES CONTRIBUTING BEYOND CDR/DXSUM\n');
    disp(AdditionalBySource);
end

%% Participant-level source evidence strings
EvidenceSources=strings(height(B),1);
for i=1:height(B)
    q=find(srcEvidence(i,:));
    if ~isempty(q)
        EvidenceSources(i)=strjoin(srcNames(q),"; ");
    end
end
B.RegistryEvidenceSources=EvidenceSources;

%% Mandatory historical classification verification
assert(height(B)==210,'Baseline cohort count mismatch.');
assert(sum(B.LongitudinalEligible)==143, ...
    'Longitudinal cohort count mismatch.');
assert(sum(nonret)==67,'Non-retained count mismatch.');
assert(nLater==22,'Broad later-activity count mismatch.');
assert(nNoLater==45,'No-later-activity count mismatch.');
assert(nCN==18,'CN later-activity count mismatch.');

fprintf('GITHUB STEP 11B VERIFICATION: PASS\n');

%% Save
writetable(B,fullfile(outputDir,'Step11B_BroadRegistry_ParticipantAudit.csv'));
writetable(SourceAudit,fullfile(outputDir,'Step11B_BroadRegistry_SourceAudit.csv'));

if ~isempty(AdditionalBySource)
    writetable(AdditionalBySource, ...
        fullfile(outputDir,'Step11B_Sources_Adding_Beyond_CDR_DXSUM.csv'));
end

% Save raw participant x source evidence for reproducibility
E=array2table(srcEvidence,'VariableNames',matlab.lang.makeValidName(cellstr(srcNames)));
E.RID=B.RID;
E.LongitudinalEligible=B.LongitudinalEligible;
E=[E(:,end-1:end) E(:,1:end-2)];
writetable(E,fullfile(outputDir,'Step11B_ParticipantBySource_RegistryEvidence.csv'));

Summary=table(height(B),sum(B.LongitudinalEligible),sum(nonret),nLater,nNoLater,nCN, ...
    'VariableNames',{'BaselineEligible','LongitudinalEligible','NonRetained', ...
    'BroadLaterRegistryActivity','NoLaterRegistryActivity', ...
    'CN_AmongBroadLaterRegistry'});
writetable(Summary,fullfile(outputDir,'Step11B_BroadRegistry_Summary.csv'));

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');

%% Local functions
function x=toNum(v)
if isnumeric(v)
    x=double(v);
else
    x=str2double(string(v));
end
end

function d=pdate(x)
if isdatetime(x)
    d=x; return
end
if isnumeric(x)
    % Numeric date values are not assumed to be MATLAB datenums because ADNI
    % fields may contain years/codes. Treat 8-digit YYYYMMDD explicitly.
    sx=string(x);
else
    sx=strtrim(string(x));
end

d=NaT(size(sx));
r=~ismissing(sx)&strlength(sx)>0;

fmts={'yyyy-MM-dd','yyyy-MM-dd HH:mm:ss','MM/dd/yyyy','M/d/yyyy', ...
      'dd-MMM-yyyy','dd/MM/yyyy','yyyyMMdd','MM/dd/yy'};

for k=1:numel(fmts)
    if ~any(r), break; end
    try
        z=datetime(sx(r),'InputFormat',fmts{k});
        ok=~isnat(z); ix=find(r);
        d(ix(ok))=z(ok); r(ix(ok))=false;
    catch
    end
end

if any(r)
    try
        z=datetime(sx(r));
        ix=find(r); ok=~isnat(z);
        d(ix(ok))=z(ok);
    catch
    end
end
d.Format='yyyy-MM-dd';
end
