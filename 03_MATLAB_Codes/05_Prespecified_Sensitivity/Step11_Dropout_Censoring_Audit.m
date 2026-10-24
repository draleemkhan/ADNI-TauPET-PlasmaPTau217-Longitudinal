%% Step11_Dropout_Censoring_Audit.m
clear; clc; close all;
sourceDir='D:\source files';
outputDir='D:\output';
F143=readtable(fullfile(outputDir,'Step08_Final143_With_Demographics.csv'),'VariableNamingRule','preserve');
CDR=readtable(fullfile(sourceDir,'CDR_11Aug2026.csv'),'VariableNamingRule','preserve');
DX=readtable(fullfile(sourceDir,'DXSUM_11Aug2026.csv'),'VariableNamingRule','preserve');

fprintf('\n====================================================\n STEP 11 - DROPOUT AND CENSORING AUDIT\n====================================================\n');

% Explicit baseline-eligible cohort dependency.
baseFile=fullfile(outputDir, ...
    'Step05A_After_AdvancedDementia_Exclusion.csv');
assert(isfile(baseFile), ...
    'Missing Step05A baseline N=210 cohort.');
B=readtable(baseFile,'VariableNamingRule','preserve');
assert(height(B)==210,'Baseline cohort must contain 210 rows.');
assert(ismember('RID',B.Properties.VariableNames), ...
    'Baseline cohort RID field missing.');
fprintf('Baseline source: %s\n',baseFile);

if ismember('IndexFTP_Date',B.Properties.VariableNames), t0=pdate(B.IndexFTP_Date);
elseif ismember('SCANDATE',B.Properties.VariableNames), t0=pdate(B.SCANDATE);
else, error('Index FTP date missing.'); end

B.LongitudinalEligible=ismember(double(B.RID),double(F143.RID));
fprintf('Baseline eligible: %d [210]\n',height(B));
fprintf('Longitudinal eligible: %d [143]\n',sum(B.LongitudinalEligible));
fprintf('No qualifying longitudinal CDR-SB: %d [67]\n',sum(~B.LongitudinalEligible));

% Registry activity: any dated CDR or DXSUM record >=3 months to <=4 years.
crid=double(CDR.RID);
if ismember('EXAMDATE',CDR.Properties.VariableNames), cd=pdate(CDR.EXAMDATE); else, cd=pdate(CDR.VISDATE); end
drid=double(DX.RID);
if ismember('EXAMDATE',DX.Properties.VariableNames), dd=pdate(DX.EXAMDATE); else, dd=pdate(DX.VISDATE); end

later=false(height(B),1);
for i=1:height(B)
    rid=double(B.RID(i));
    a=days(cd(crid==rid)-t0(i))/365.25;
    b=days(dd(drid==rid)-t0(i))/365.25;
    later(i)=any(a>=.25 & a<=4) || any(b>=.25 & b<=4);
end
B.LaterRegistryActivity=later;
nr=~B.LongitudinalEligible;
fprintf('\nLater registry activity, no qualifying CDR-SB: %d [verified 21]\n',sum(nr&later));
fprintf('No later registry activity: %d [verified 46]\n',sum(nr&~later));

% Baseline diagnosis nearest +/-90 d.
raw=DX.DIAGNOSIS;
if isnumeric(raw), dc=double(raw); else, dc=str2double(string(raw)); end
bdx=nan(height(B),1);
for i=1:height(B)
    q=find(drid==double(B.RID(i)) & ~isnat(dd) & isfinite(dc));
    if isempty(q), continue; end
    g=days(dd(q)-t0(i)); [m,j]=min(abs(g));
    if m<=90, bdx(i)=dc(q(j)); end
end
B.BaselineDiagnosisCode=bdx;
B.BaselineDiagnosis=repmat("Unclassified",height(B),1);
B.BaselineDiagnosis(bdx==1)="Cognitively normal";
B.BaselineDiagnosis(bdx==2)="MCI";
B.BaselineDiagnosis(bdx==3)="Dementia";
fprintf('CN among later-registry/non-retained: %d / %d [verified 17/21]\n', ...
    sum(nr&later&B.BaselineDiagnosis=="Cognitively normal"),sum(nr&later));

% Follow-up opportunity: use latest dated CDR/DXSUM record as administrative cutoff.
allD=[cd;dd]; allD=allD(~isnat(allD)); cutoff=max(allD);
B.FollowupOpportunityYears=days(cutoff-t0)/365.25;
B.Adequate4YearOpportunity=B.FollowupOpportunityYears>=4;
fprintf('\nAdministrative cutoff from source records: %s\n',string(cutoff));
fprintf('Median opportunity, retained: %.2f y\n',median(B.FollowupOpportunityYears(B.LongitudinalEligible),'omitnan'));
fprintf('Median opportunity, non-retained: %.2f y\n',median(B.FollowupOpportunityYears(nr),'omitnan'));
fprintf('Participants with >=4-year opportunity: %d / %d\n',sum(B.Adequate4YearOpportunity),height(B));

%% Mandatory verification against locked historical Step 11 summary
assert(height(B)==210,'Baseline count mismatch.');
assert(sum(B.LongitudinalEligible)==143, ...
    'Longitudinal count mismatch.');
assert(sum(nr)==67,'Non-retained count mismatch.');
assert(sum(nr&later)==21,'Later registry count mismatch.');
assert(sum(nr&~later)==46,'No later registry count mismatch.');
assert(sum(nr&later&B.BaselineDiagnosis=="Cognitively normal")==17, ...
    'CN among later registry participants mismatch.');
fprintf('GITHUB STEP 11 VERIFICATION: PASS\n');

% Save participant-level audit first; detailed baseline comparisons follow once
% we confirm exact variable names in the N=210 source file.
writetable(B,fullfile(outputDir,'Step11_Baseline210_Dropout_Censoring_Audit.csv'));

S=table(height(B),sum(B.LongitudinalEligible),sum(nr),sum(nr&later),sum(nr&~later), ...
    sum(nr&later&B.BaselineDiagnosis=="Cognitively normal"), ...
    'VariableNames',{'BaselineEligible','LongitudinalEligible','NoQualifyingCDRSB', ...
    'LaterRegistry_NoQualifyingCDRSB','NoLaterRegistryActivity','CN_LaterRegistryNonRetained'});
writetable(S,fullfile(outputDir,'Step11_Dropout_Censoring_Summary.csv'));

fprintf('\nAvailable variables in selected N=210 baseline file:\n');
disp(string(B.Properties.VariableNames)');
fprintf('\nOutputs saved to:\n%s\n====================================================\n',outputDir);

function d=pdate(x)
if isdatetime(x), d=x; return; end
s=strtrim(string(x)); d=NaT(size(s));
f={'yyyy-MM-dd','yyyy-MM-dd HH:mm:ss','MM/dd/yyyy','MM/yyyy','dd-MMM-yyyy','dd/MM/yyyy'};
r=~ismissing(s)&strlength(s)>0;
for k=1:numel(f)
    if ~any(r), break; end
    try
        z=datetime(s(r),'InputFormat',f{k}); ok=~isnat(z); ix=find(r);
        d(ix(ok))=z(ok); r(ix(ok))=false;
    catch
    end
end
if any(r)
    try
        z=datetime(s(r)); ix=find(r); ok=~isnat(z); d(ix(ok))=z(ok);
    catch
    end
end
d.Format='yyyy-MM-dd';
end
