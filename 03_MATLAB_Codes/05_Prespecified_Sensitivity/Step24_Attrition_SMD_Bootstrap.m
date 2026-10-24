%% Step24_Attrition_SMD_Bootstrap.m
% ADNI revised cohort model
% Step 24: Attrition/retention comparison using standardized mean
% differences (SMDs) with participant-level bootstrap 95% CIs.
%
% Input:
% D:\output\
% Step11C_Baseline210_With_Demographics.csv
%
% Analyses:
%   A) Full baseline-eligible cohort: retained vs non-retained
%   B) Participants with >=4 years follow-up opportunity
%
% Continuous variables: standardized mean difference using pooled SD.
% Binary variables: difference in proportions standardized by pooled
% Bernoulli SD.
% Diagnosis is represented by separate binary indicators (CN/MCI/Dementia).
%
% Bootstrap: 2000 stratified resamples WITHIN retained/non-retained groups.
% This preserves the observed group sizes and estimates percentile 95% CIs.
%
% IMPORTANT:
% SMD direction = Retained minus Non-retained.
% Positive SMD therefore means higher/more frequent in retained participants.

clear; clc;

%% Paths/settings
outDir = 'D:\output';
B = 2000;
rng(24092026,'twister');

% Accept either the original local filename or the "(1)" uploaded-copy name.
candidates = { ...
    fullfile(outDir,'Step11C_Baseline210_With_Demographics.csv'), ...
    fullfile(outDir,'Step11C_Baseline210_With_Demographics(1).csv')};

inFile = '';
for k = 1:numel(candidates)
    if isfile(candidates{k})
        inFile = candidates{k};
        break;
    end
end
assert(~isempty(inFile), ...
    'Could not find Step11C_Baseline210_With_Demographics.csv in D:\output.');

fprintf('\nSTEP 24: ATTRITION / RETENTION SMD ANALYSIS\n');
fprintf('===========================================\n');
fprintf('Input: %s\n',inFile);
fprintf('Bootstrap replicates: %d\n',B);

T = readtable(inFile,'VariableNamingRule','preserve');
fprintf('Loaded: %d participants, %d variables\n',height(T),width(T));

%% Resolve variable names robustly
vRet = findVar(T, {'LongitudinalEligible'});
vOpp = findVar(T, {'FollowupOpportunityYears','FollowUpOpportunityYears'});
vCDR = findVarContains(T, {'CDRSB','CDR_SB'}, {'baseline','base'}, true);
vTau = findVarContains(T, {'TauMetaSUVR','FTP','MetaTemporal','Meta_Temporal'}, {}, true);
vCL  = findVarContains(T, {'Centiloid'}, {}, true);
vPT  = findVarContains(T, {'PTau217','pTau217','pT217'}, {}, true);
vAge = findVar(T, {'Age'});
vEdu = findVar(T, {'Education','PTEDUCAT'});
vSex = findVar(T, {'Sex','PTGENDER'});
vDx  = findVarContains(T, {'Diagnosis','DX'}, {'baseline','base'}, false);

fprintf('\nResolved variables:\n');
fprintf('Retention: %s\n',vRet);
fprintf('Opportunity: %s\n',vOpp);
fprintf('CDR-SB: %s\n',vCDR);
fprintf('Tau PET: %s\n',vTau);
fprintf('Centiloid: %s\n',vCL);
fprintf('p-tau217: %s\n',vPT);
fprintf('Age: %s\n',vAge);
fprintf('Education: %s\n',vEdu);
fprintf('Sex: %s\n',vSex);
fprintf('Diagnosis: %s\n',vDx);

%% Convert retention flag
ret = toBinaryRetention(T.(vRet));
assert(sum(ret==1)==143 && sum(ret==0)==67, ...
    'Expected 143 retained and 67 non-retained in full N=210 cohort.');

%% Create clean analysis variables
A = table;
A.Retained = ret;
A.CDRSB = toNumeric(T.(vCDR));
A.TauPET = toNumeric(T.(vTau));
A.Centiloid = toNumeric(T.(vCL));
A.PTau217 = toNumeric(T.(vPT));
A.Age = toNumeric(T.(vAge));
A.Education = toNumeric(T.(vEdu));
A.Female = sexToFemale(T.(vSex));
A.FollowupOpportunityYears = toNumeric(T.(vOpp));

[cn,mci,dem] = diagnosisIndicators(T.(vDx));
A.CN = cn;
A.MCI = mci;
A.Dementia = dem;

%% Audit
fprintf('\nFULL COHORT AUDIT\n');
fprintf('N = %d | retained = %d | non-retained = %d\n', ...
    height(A),sum(A.Retained==1),sum(A.Retained==0));
fprintf('Missing values by analysis variable:\n');
disp(varfun(@(x)sum(isnan(x)),A,'OutputFormat','table'));

%% Full cohort SMDs
fprintf('\nA) FULL BASELINE-ELIGIBLE COHORT\n');
Full210 = buildSMDTable(A,B);
disp(Full210);

%% >=4-year opportunity subset
A4 = A(A.FollowupOpportunityYears >= 4,:);
fprintf('\nB) >=4-YEAR FOLLOW-UP-OPPORTUNITY SUBSET\n');
fprintf('N = %d | retained = %d | non-retained = %d\n', ...
    height(A4),sum(A4.Retained==1),sum(A4.Retained==0));

assert(height(A4)==125 && sum(A4.Retained==1)==94 && sum(A4.Retained==0)==31, ...
    ['The >=4-year subset does not reproduce expected N=125 ' ...
     '(94 retained, 31 non-retained). Check FollowupOpportunityYears.']);

Opportunity4Y = buildSMDTable(A4,B);
disp(Opportunity4Y);

%% Compact combined manuscript-oriented table
FullCompact = Full210(:,{'Variable','RetainedSummary','NonRetainedSummary', ...
    'SMD','SMD_CI_Lower','SMD_CI_Upper','AbsSMD'});
FullCompact.Analysis = repmat("Full baseline-eligible cohort",height(FullCompact),1);

FourCompact = Opportunity4Y(:,{'Variable','RetainedSummary','NonRetainedSummary', ...
    'SMD','SMD_CI_Lower','SMD_CI_Upper','AbsSMD'});
FourCompact.Analysis = repmat(">=4-year follow-up opportunity",height(FourCompact),1);

Combined = [FullCompact;FourCompact];
Combined = movevars(Combined,'Analysis','Before','Variable');

%% Descriptive flags (NOT hypothesis-test declarations)
% Conventional descriptive thresholds only; interpretation should focus
% on magnitude and CIs rather than treating 0.1 as a significance cutoff.
Combined.AbsSMD_ge_0_10 = Combined.AbsSMD >= 0.10;
Combined.AbsSMD_ge_0_20 = Combined.AbsSMD >= 0.20;

fprintf('\nLARGEST ABSOLUTE SMDs\n');
tmp1 = sortrows(Full210,'AbsSMD','descend');
fprintf('\nFull N=210:\n');
disp(tmp1(:,{'Variable','SMD','SMD_CI_Lower','SMD_CI_Upper','AbsSMD'}));

tmp2 = sortrows(Opportunity4Y,'AbsSMD','descend');
fprintf('\n>=4-year opportunity subset:\n');
disp(tmp2(:,{'Variable','SMD','SMD_CI_Lower','SMD_CI_Upper','AbsSMD'}));

%% Save
matFile = fullfile(outDir,'Step24_Attrition_SMD_Bootstrap.mat');
save(matFile,'A','A4','Full210','Opportunity4Y','Combined','B','inFile','-v7.3');

writetable(Full210,fullfile(outDir,'Step24_Attrition_SMD_Full210.csv'));
writetable(Opportunity4Y,fullfile(outDir,'Step24_Attrition_SMD_4YearOpportunity.csv'));
writetable(Combined,fullfile(outDir,'Step24_Attrition_SMD_Combined.csv'));

fprintf('\nSTEP 24 COMPLETE\n');
fprintf('Saved outputs to:\n%s\n',outDir);
fprintf('\nPlease paste:\n');
fprintf('1) FULL COHORT table,\n');
fprintf('2) >=4-YEAR table, and\n');
fprintf('3) LARGEST ABSOLUTE SMD tables.\n');

%% ================= LOCAL FUNCTIONS =================

function R = buildSMDTable(A,B)
    vars = {'CDRSB','TauPET','Centiloid','PTau217','Age','Education', ...
            'Female','CN','MCI','Dementia'};
    labels = ["Baseline CDR-SB","Tau PET meta-temporal SUVR","Centiloid", ...
              "Plasma p-tau217","Age","Education","Female sex", ...
              "CN diagnosis","MCI diagnosis","Dementia diagnosis"];
    isBinary = [false false false false false false true true true true];

    nV = numel(vars);
    SMD = nan(nV,1); lo = nan(nV,1); hi = nan(nV,1);
    AbsSMD = nan(nV,1);
    RetainedSummary = strings(nV,1);
    NonRetainedSummary = strings(nV,1);
    NRetained = nan(nV,1); NNonRetained = nan(nV,1);

    g1idx = find(A.Retained==1);
    g0idx = find(A.Retained==0);

    for j = 1:nV
        x = A.(vars{j});
        x1 = x(g1idx); x1 = x1(isfinite(x1));
        x0 = x(g0idx); x0 = x0(isfinite(x0));
        NRetained(j)=numel(x1); NNonRetained(j)=numel(x0);

        if isBinary(j)
            SMD(j)=smdBinary(x1,x0);
            RetainedSummary(j)=sprintf('%d/%d (%.1f%%)',sum(x1==1),numel(x1),100*mean(x1));
            NonRetainedSummary(j)=sprintf('%d/%d (%.1f%%)',sum(x0==1),numel(x0),100*mean(x0));
        else
            SMD(j)=smdContinuous(x1,x0);
            RetainedSummary(j)=sprintf('%.3f (SD %.3f)',mean(x1),std(x1));
            NonRetainedSummary(j)=sprintf('%.3f (SD %.3f)',mean(x0),std(x0));
        end

        % Stratified bootstrap within retention groups.
        bs = nan(B,1);
        for b=1:B
            b1=x1(randi(numel(x1),numel(x1),1));
            b0=x0(randi(numel(x0),numel(x0),1));
            if isBinary(j)
                bs(b)=smdBinary(b1,b0);
            else
                bs(b)=smdContinuous(b1,b0);
            end
        end
        bs=bs(isfinite(bs));
        if ~isempty(bs)
            ci=prctile(bs,[2.5 97.5]);
            lo(j)=ci(1); hi(j)=ci(2);
        end
        AbsSMD(j)=abs(SMD(j));
    end

    R=table(labels',NRetained,NNonRetained,RetainedSummary,NonRetainedSummary, ...
        SMD,lo,hi,AbsSMD, ...
        'VariableNames',{'Variable','N_Retained','N_NonRetained', ...
        'RetainedSummary','NonRetainedSummary','SMD','SMD_CI_Lower', ...
        'SMD_CI_Upper','AbsSMD'});
end

function d=smdContinuous(x1,x0)
    n1=numel(x1); n0=numel(x0);
    if n1<2 || n0<2, d=NaN; return; end
    sp=sqrt(((n1-1)*var(x1)+(n0-1)*var(x0))/(n1+n0-2));
    if sp==0, d=0; else, d=(mean(x1)-mean(x0))/sp; end
end

function d=smdBinary(x1,x0)
    p1=mean(x1); p0=mean(x0);
    p=(p1+p0)/2;
    den=sqrt(p*(1-p));
    if den==0, d=0; else, d=(p1-p0)/den; end
end

function name=findVar(T,candidates)
    vn=T.Properties.VariableNames;
    name='';
    for i=1:numel(candidates)
        k=find(strcmpi(vn,candidates{i}),1);
        if ~isempty(k), name=vn{k}; return; end
    end
    error('Could not locate variable. Candidates: %s',strjoin(candidates,', '));
end

function name=findVarContains(T,keys,preferred,requirePreferred)
    vn=T.Properties.VariableNames;
    score=zeros(size(vn));
    for i=1:numel(vn)
        s=lower(vn{i});
        for k=1:numel(keys)
            if contains(s,lower(keys{k})), score(i)=score(i)+10; end
        end
        for k=1:numel(preferred)
            if contains(s,lower(preferred{k})), score(i)=score(i)+2; end
        end
    end
    [mx,idx]=max(score);
    if mx<10
        error('Could not locate variable containing: %s',strjoin(keys,', '));
    end
    if requirePreferred && ~isempty(preferred)
        % Prefer baseline-labelled match if one exists, but do not require it
        % when only one unambiguous biomarker variable is present.
    end
    name=vn{idx};
end

function x=toNumeric(v)
    if isnumeric(v)
        x=double(v);
    elseif islogical(v)
        x=double(v);
    elseif iscategorical(v)
        x=str2double(string(v));
    else
        x=str2double(string(v));
    end
end

function r=toBinaryRetention(v)
    if isnumeric(v) || islogical(v)
        r=double(v);
        return;
    end
    s=lower(strtrim(string(v)));
    r=nan(numel(s),1);
    yes=ismember(s,["1","true","yes","y","eligible","retained"]);
    no=ismember(s,["0","false","no","n","ineligible","non-retained","nonretained"]);
    r(yes)=1; r(no)=0;
    assert(all(isfinite(r)),'Could not convert LongitudinalEligible to binary.');
end

function f=sexToFemale(v)
    if isnumeric(v) || islogical(v)
        u=unique(v(~isnan(v)));
        if all(ismember(u,[0 1]))
            % Cannot infer coding safely.
            error('Sex is numeric 0/1; coding cannot be inferred safely.');
        end
    end
    s=upper(strtrim(string(v)));
    f=nan(numel(s),1);
    f(ismember(s,["F","FEMALE","2"]))=1;
    f(ismember(s,["M","MALE","1"]))=0;
    assert(all(isfinite(f)),'Could not convert sex variable to Female indicator.');
end

function [cn,mci,dem]=diagnosisIndicators(v)
    n=numel(v); cn=nan(n,1); mci=nan(n,1); dem=nan(n,1);

    if isnumeric(v)
        x=double(v);
        cn=double(x==1); mci=double(x==2); dem=double(x==3);
        bad=~ismember(x,[1 2 3]) & ~isnan(x);
        assert(~any(bad),'Unexpected numeric diagnosis codes.');
        cn(isnan(x))=NaN; mci(isnan(x))=NaN; dem(isnan(x))=NaN;
        return;
    end

    s=upper(strtrim(string(v)));
    cn(:)=0; mci(:)=0; dem(:)=0;
    cn(contains(s,"CN") | contains(s,"NORMAL"))=1;
    mci(contains(s,"MCI"))=1;
    dem(contains(s,"DEM") | contains(s,"AD"))=1;
    known=(cn+mci+dem)>0;
    cn(~known)=NaN; mci(~known)=NaN; dem(~known)=NaN;
    assert(all(known),'Could not classify one or more baseline diagnoses.');
end
