%% Step11C_Retained_vs_NonRetained_Attrition_Comparison_v3.m
% ADNI Tau PET revised cohort analysis
% Fixes MATLAB indexed table-variable syntax in tau-specific attrition check.
% Statistical comparisons are unchanged from v2.

clear; clc; close all;

sourceDir = 'D:\source files';
outputDir = 'D:\output';

auditFile = fullfile(outputDir,'Step11B_BroadRegistry_ParticipantAudit.csv');
demogFile = fullfile(sourceDir,'PTDEMOG_11Aug2026.csv');
cdrFile   = fullfile(sourceDir,'CDR_11Aug2026.csv');
dxFile    = fullfile(sourceDir,'DXSUM_11Aug2026.csv');

assert(isfile(auditFile),'Missing Step11B participant audit.');
assert(isfile(demogFile),'Missing PTDEMOG source file.');
assert(isfile(cdrFile),'Missing CDR source file.');
assert(isfile(dxFile),'Missing DXSUM source file.');

B = readtable(auditFile,'VariableNamingRule','preserve');
DEM = readtable(demogFile,'VariableNamingRule','preserve');
CDR = readtable(cdrFile,'VariableNamingRule','preserve');
DX  = readtable(dxFile,'VariableNamingRule','preserve');

fprintf('\n====================================================\n');
fprintf(' STEP 11C v3 - RETAINED VS NON-RETAINED COMPARISON\n');
fprintf('====================================================\n');

%% Index date
if ismember('IndexFTP_Date_Audit',B.Properties.VariableNames)
    t0 = pdate(B.IndexFTP_Date_Audit);
elseif ismember('IndexFTP_Date',B.Properties.VariableNames)
    t0 = pdate(B.IndexFTP_Date);
elseif ismember('SCANDATE',B.Properties.VariableNames)
    t0 = pdate(B.SCANDATE);
else
    error('Index FTP date missing.');
end

%% Follow-up opportunity
if ~ismember('FollowupOpportunityYears',B.Properties.VariableNames)
    if ismember('EXAMDATE',CDR.Properties.VariableNames)
        cdrDate = pdate(CDR.EXAMDATE);
    else
        cdrDate = pdate(CDR.VISDATE);
    end

    if ismember('EXAMDATE',DX.Properties.VariableNames)
        dxDate = pdate(DX.EXAMDATE);
    else
        dxDate = pdate(DX.VISDATE);
    end

    allDates=[cdrDate;dxDate];
    allDates=allDates(~isnat(allDates));
    cutoff=max(allDates);
    B.FollowupOpportunityYears=days(cutoff-t0)/365.25;
else
    cutoff=NaT;
end

%% Demographics
dRID=toNum(DEM.RID);
if ismember('VISDATE',DEM.Properties.VariableNames)
    dDate=pdate(DEM.VISDATE);
else
    dDate=NaT(height(DEM),1);
end
dDOB=pdate(DEM.PTDOBYY);
dSex=toNum(DEM.PTGENDER);
dEdu=toNum(DEM.PTEDUCAT);

Age=nan(height(B),1);
Sex=repmat("Missing",height(B),1);
Education=nan(height(B),1);

for i=1:height(B)
    rid=double(B.RID(i));
    q=find(dRID==rid);

    qAS=q(~isnat(dDOB(q)) & ismember(dSex(q),[1 2]));
    if ~isempty(qAS)
        qDated=qAS(~isnat(dDate(qAS)));
        if ~isempty(qDated)
            [~,j]=min(dDate(qDated)); k=qDated(j);
        else
            k=qAS(1);
        end
        Age(i)=years(t0(i)-dDOB(k));
        if dSex(k)==1, Sex(i)="Male"; else, Sex(i)="Female"; end
    end

    qE=q(isfinite(dEdu(q)) & dEdu(q)>=0 & dEdu(q)<=30);
    if ~isempty(qE)
        qEDated=qE(~isnat(dDate(qE)));
        if ~isempty(qEDated)
            [~,j]=min(dDate(qEDated)); k=qEDated(j);
        else
            k=qE(1);
        end
        Education(i)=dEdu(k);
    end
end

B.Age=Age;
B.Sex=Sex;
B.Education=Education;

tauVar='META_TEMPORAL_SUVR';
centVar='Centiloid';
ptauVar='PTau217_Fujirebio';
cdrVar='Baseline_CDRSB';

%% Full comparison
fullMask=true(height(B),1);
Full=makeComparison(B,fullMask,tauVar,centVar,ptauVar,cdrVar);

fprintf('\nFULL BASELINE-ELIGIBLE COHORT (N=%d)\n',height(B));
disp(Full);

%% Adequate opportunity comparison
adequateMask=B.FollowupOpportunityYears>=4;
Adequate=makeComparison(B,adequateMask,tauVar,centVar,ptauVar,cdrVar);

fprintf('\nRESTRICTED TO >=4 YEARS FOLLOW-UP OPPORTUNITY\n');
fprintf('Participants included: %d\n',sum(adequateMask));
fprintf('Retained: %d | Non-retained: %d\n', ...
    sum(adequateMask & B.LongitudinalEligible), ...
    sum(adequateMask & ~B.LongitudinalEligible));
disp(Adequate);

%% Tau-specific attrition check - fixed syntax
tauAll = double(B.(tauVar));
retFlag = logical(B.LongitudinalEligible);

tauRet = tauAll(retFlag & isfinite(tauAll));
tauNo  = tauAll(~retFlag & isfinite(tauAll));
pTauFull = ranksum(tauRet,tauNo);

tauRetA = tauAll(adequateMask & retFlag & isfinite(tauAll));
tauNoA  = tauAll(adequateMask & ~retFlag & isfinite(tauAll));

if ~isempty(tauRetA) && ~isempty(tauNoA)
    pTauAdeq = ranksum(tauRetA,tauNoA);
else
    pTauAdeq = NaN;
end

fprintf('\nTAU-SPECIFIC ATTRITION CHECK\n');
fprintf('Full cohort Tau difference P = %.6g\n',pTauFull);
fprintf('>=4-year opportunity Tau difference P = %.6g\n',pTauAdeq);

fprintf('\nFOLLOW-UP OPPORTUNITY\n');
fprintf('Median retained: %.2f y\n',median(B.FollowupOpportunityYears(retFlag),'omitnan'));
fprintf('Median non-retained: %.2f y\n',median(B.FollowupOpportunityYears(~retFlag),'omitnan'));

%% Mandatory verification against historical Step 11C summary
assert(height(B)==210,'Baseline cohort mismatch.');
assert(sum(retFlag)==143,'Retained count mismatch.');
assert(sum(~retFlag)==67,'Non-retained count mismatch.');
assert(sum(adequateMask)==125, ...
    'Four-year opportunity count mismatch.');
assert(sum(adequateMask&retFlag)==94, ...
    'Four-year retained count mismatch.');
assert(sum(adequateMask&~retFlag)==31, ...
    'Four-year non-retained count mismatch.');
assert(isfinite(pTauFull) && abs(pTauFull-0.15478)<1e-5, ...
    'Full-cohort tau P-value mismatch.');
assert(isfinite(pTauAdeq) && abs(pTauAdeq-0.67017)<1e-5, ...
    'Four-year tau P-value mismatch.');
fprintf('GITHUB STEP 11C VERIFICATION: PASS\n');

%% Save outputs
writetable(B,fullfile(outputDir,'Step11C_Baseline210_With_Demographics.csv'));
writetable(Full,fullfile(outputDir,'Step11C_AttritionComparison_Full210.csv'));
writetable(Adequate,fullfile(outputDir,'Step11C_AttritionComparison_4YearOpportunity.csv'));

Summary=table( ...
    height(B),sum(retFlag),sum(~retFlag),sum(adequateMask), ...
    sum(adequateMask&retFlag),sum(adequateMask&~retFlag), ...
    pTauFull,pTauAdeq, ...
    'VariableNames',{'BaselineEligible','Retained','NonRetained', ...
    'AdequateOpportunityN','AdequateRetained','AdequateNonRetained', ...
    'TauP_Full','TauP_AdequateOpportunity'});
writetable(Summary,fullfile(outputDir,'Step11C_AttritionComparison_Summary.csv'));

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');

%% Local functions
function Out=makeComparison(T,mask,tauVar,centVar,ptauVar,cdrVar)
R=logical(T.LongitudinalEligible)&mask;
N=~logical(T.LongitudinalEligible)&mask;

labels=[ ...
"Baseline CDR-SB"
"Tau PET meta-temporal SUVR"
"Amyloid burden, Centiloids"
"Plasma p-tau217"
"Age, years"
"Education, years"
"Female sex, n (%)"
"Cognitively normal, n (%)"
"MCI, n (%)"
"Dementia, n (%)"];

Retained=strings(numel(labels),1);
NonRetained=strings(numel(labels),1);
PValue=nan(numel(labels),1);

vars={cdrVar,tauVar,centVar,ptauVar,'Age','Education'};

for i=1:6
    x=double(T.(vars{i}));
    xR=x(R&isfinite(x));
    xN=x(N&isfinite(x));
    Retained(i)=fmtMed(xR);
    NonRetained(i)=fmtMed(xN);
    if ~isempty(xR)&&~isempty(xN), PValue(i)=ranksum(xR,xN); end
end

flags={string(T.Sex)=="Female", ...
       string(T.BaselineDiagnosis)=="Cognitively normal", ...
       string(T.BaselineDiagnosis)=="MCI", ...
       string(T.BaselineDiagnosis)=="Dementia"};

for j=1:4
    i=6+j; f=flags{j};
    Retained(i)=fmtPct(f(R),sum(R));
    NonRetained(i)=fmtPct(f(N),sum(N));
    tab=[sum(R&f),sum(R&~f);sum(N&f),sum(N&~f)];
    try
        [~,PValue(i)]=fishertest(tab);
    catch
        PValue(i)=NaN;
    end
end

Out=table(labels,Retained,NonRetained,PValue, ...
    'VariableNames',{'Characteristic','Retained','NonRetained','PValue'});
end

function s=fmtMed(x)
if isempty(x), s="NA"; return; end
q=prctile(x,[25 75]);
s=sprintf('%.3f (%.3f-%.3f)',median(x),q(1),q(2));
end

function s=fmtPct(flag,N)
n=sum(flag);
if N==0, s="NA"; else, s=sprintf('%d (%.1f%%)',n,100*n/N); end
end

function x=toNum(v)
if isnumeric(v), x=double(v); else, x=str2double(string(v)); end
end

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
