%% Step16_StageSpecific_Tau_Sensitivity.m
% ADNI Tau PET revised cohort analysis
%
% STEP 16 — Stage-specific tau sensitivity
%
% Fits the LOCKED uncorrelated random-intercept + random-Time-slope LMM
% separately within baseline diagnostic groups:
%   Cognitively normal, MCI, Dementia
%
% Expected:
%   CN:       33 participants, 84 observations
%   MCI:      75 participants, 208 observations
%   Dementia: 35 participants, 87 observations
%
% Historical Tau x Time benchmarks:
%   CN        beta ~0.017, P ~.905
%   MCI       beta ~0.172, P ~.136
%   Dementia  beta ~1.105, P ~.022
%
% IMPORTANT:
% These are exploratory subgroup sensitivity analyses. Differences in
% subgroup P values do NOT establish statistically significant effect
% modification by diagnosis.

clear; clc; close all;

outputDir = 'D:\output';
dataFile  = fullfile(outputDir,'Step08_ModelReady_379_Observations.csv');
dxFile    = fullfile(outputDir,'Step07_Final143_With_Baseline_Diagnosis.csv');

assert(isfile(dataFile),'Missing Step08 model-ready longitudinal dataset.');
assert(isfile(dxFile),'Missing Step07 diagnosis dataset.');

T = readtable(dataFile,'VariableNamingRule','preserve');
D = readtable(dxFile,'VariableNamingRule','preserve');

fprintf('\n====================================================\n');
fprintf(' STEP 16 - STAGE-SPECIFIC TAU SENSITIVITY\n');
fprintf('====================================================\n');

%% Ensure numeric variables
numVars={'RID','CDRSB','TimeYears','Tau_z','Centiloid_z','PTau217_z','Age','Education'};
for k=1:numel(numVars)
    if ~isnumeric(T.(numVars{k}))
        T.(numVars{k})=str2double(string(T.(numVars{k})));
    end
end

%% Attach diagnosis
if ~ismember('Baseline_Diagnosis',T.Properties.VariableNames)
    if ~all(ismember({'RID','Baseline_Diagnosis'},D.Properties.VariableNames))
        error('Baseline diagnosis unavailable.');
    end
    T=innerjoin(T,D(:,{'RID','Baseline_Diagnosis'}),'Keys','RID');
end

T.RIDcat=categorical(T.RID);
T.SexModel=categorical(string(T.Sex),["Male","Female"]);

complete=isfinite(T.CDRSB) & isfinite(T.TimeYears) & ...
         isfinite(T.Tau_z) & isfinite(T.Centiloid_z) & ...
         isfinite(T.PTau217_z) & isfinite(T.Age) & ...
         isfinite(T.Education) & ~isundefined(T.SexModel);

T=T(complete,:);

fprintf('Overall complete observations: %d [379]\n',height(T));
fprintf('Overall participants: %d [143]\n',numel(unique(T.RID)));

%% Same fixed-effects specification as primary model
formula=['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
         'Age + SexModel + Education + Tau_z:TimeYears + ' ...
         'Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
         '(1 + TimeYears|RIDcat)'];

groups=["Cognitively normal","MCI","Dementia"];

Result=table();
AllFE=cell(3,1);

for g=1:numel(groups)
    S=T(string(T.Baseline_Diagnosis)==groups(g),:);

    nObs=height(S);
    nSub=numel(unique(S.RID));

    fprintf('\n----------------------------------------------------\n');
    fprintf('%s\n',upper(groups(g)));
    fprintf('Participants: %d\n',nSub);
    fprintf('CDR-SB observations: %d\n',nObs);
    fprintf('----------------------------------------------------\n');

    try
        mdl=fitlme(S,formula,'FitMethod','REML', ...
            'CovariancePattern','Diagonal');

        FEraw=mdl.Coefficients;
        if istable(FEraw)
            FE=FEraw;
        else
            FE=dataset2table(FEraw);
        end
        AllFE{g}=FE;

        nm=string(FE.Name);
        ix=find(nm=="TimeYears:Tau_z" | nm=="Tau_z:TimeYears",1);
        if isempty(ix), error('Tau x Time term not found.'); end

        beta=FE.Estimate(ix);
        se=FE.SE(ix);
        lo=FE.Lower(ix);
        hi=FE.Upper(ix);
        p=FE.pValue(ix);

        fprintf('Tau x Time beta = %.6f\n',beta);
        fprintf('SE = %.6f\n',se);
        fprintf('95%% CI = %.6f to %.6f\n',lo,hi);
        fprintf('P = %.8g\n',p);

        row=table(groups(g),nSub,nObs,beta,se,lo,hi,p, ...
            mdl.ModelCriterion.AIC,mdl.ModelCriterion.BIC, ...
            'VariableNames',{'BaselineDiagnosis','Participants','Observations', ...
            'TauTime_Beta','SE','CI_Lower','CI_Upper','PValue','AIC','BIC'});
    catch ME
        warning('Model failed for %s: %s',groups(g),ME.message);
        row=table(groups(g),nSub,nObs,NaN,NaN,NaN,NaN,NaN,NaN,NaN, ...
            'VariableNames',{'BaselineDiagnosis','Participants','Observations', ...
            'TauTime_Beta','SE','CI_Lower','CI_Upper','PValue','AIC','BIC'});
    end

    Result=[Result;row]; %#ok<AGROW>
end

fprintf('\n====================================================\n');
fprintf(' STAGE-SPECIFIC TAU SENSITIVITY SUMMARY\n');
fprintf('====================================================\n');
disp(Result);

%% Save
writetable(Result,fullfile(outputDir,'Step16_StageSpecific_Tau_Sensitivity.csv'));

for g=1:3
    if ~isempty(AllFE{g})
        safeName=regexprep(groups(g),'[^A-Za-z0-9]','_');
        writetable(AllFE{g},fullfile(outputDir, ...
            ['Step16_' char(safeName) '_FixedEffects.csv']));
    end
end

fprintf('\nInterpret subgroup estimates cautiously: subgroup-specific P values\n');
fprintf('do not constitute a formal test of effect modification by diagnosis.\n');
fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');
