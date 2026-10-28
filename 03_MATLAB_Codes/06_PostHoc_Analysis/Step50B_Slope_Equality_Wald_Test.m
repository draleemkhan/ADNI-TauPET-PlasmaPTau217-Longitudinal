%% Step50B_Slope_Equality_Wald_Test.m
% ADNI revised cohort - post hoc slope-equality test (MATLAB R2025a)
%
% PURPOSE
% Post hoc Wald test of equality of the mutually adjusted tau PET x time
% and plasma p-tau217 x time coefficients, added during manuscript
% revision. H0: beta(Tau_z:TimeYears) - beta(PTau217_z:TimeYears) = 0.
%
% The difference, its standard error, and its 95% CI use the full
% fixed-effect covariance matrix of the same fitted model, including the
% covariance between the two coefficients. The P value is from coefTest
% (F test with residual degrees of freedom), the same reference
% distribution used for fitlme coefficient P values elsewhere.
%
% MODELS
%   A. Whole cohort, primary mutually adjusted model, REML (reported)
%   B. Whole cohort, same model, ML (consistency check)
%   C. MCI subgroup, mutually adjusted model as in Steps 47/49I, ML
% Random effects: uncorrelated random intercept and random linear-time
% slope, specified as separate terms.
%
% BENCHMARKS (manuscript, Table S3)
%   A. Whole cohort REML: tau 0.413, p-tau217 0.211;
%      difference 0.202 (95% CI -0.413 to 0.817); F(1,368)=0.418; P=.52
%   B. Whole cohort ML:   difference 0.200 (95% CI -0.406 to 0.806); P=.52
%   C. MCI ML:            tau 0.172, p-tau217 0.148;
%      difference 0.023 (95% CI -0.405 to 0.451); F(1,197)=0.012; P=.92
%
% INPUT : D:\output\Step08_ModelReady_379_Observations.csv
% OUTPUT: D:\output\Step50B_Slope_Equality_Wald_Test.csv (summary only)
%
% No participant-level data are written by this script.

clear; clc;

fprintf('\n============================================================\n');
fprintf(' STEP 50B: POST HOC SLOPE-EQUALITY WALD TEST\n');
fprintf('============================================================\n');

%% Paths
outputDir = 'D:\output';
dataFile  = fullfile(outputDir,'Step08_ModelReady_379_Observations.csv');
outFile   = fullfile(outputDir,'Step50B_Slope_Equality_Wald_Test.csv');
assert(isfile(dataFile),'Step08 model-ready dataset not found.');

%% Load and prepare
T = readtable(dataFile,'VariableNamingRule','preserve');

numVars = {'RID','CDRSB','TimeYears','Tau_z','Centiloid_z','PTau217_z','Age','Education'};
for k = 1:numel(numVars)
    if ~isnumeric(T.(numVars{k}))
        T.(numVars{k}) = str2double(string(T.(numVars{k})));
    end
end
T.RIDcat   = categorical(T.RID);
T.SexModel = categorical(string(T.Sex),["Male","Female"]);   % male reference

if ismember('Baseline_Diagnosis',T.Properties.VariableNames)
    dx = string(T.Baseline_Diagnosis);
else
    dx = string(T.Diagnosis);
end

assert(height(T)==379 && numel(unique(T.RID))==143, ...
    'Expected 379 observations from 143 participants.');
TMCI = T(dx=="MCI",:);
assert(height(TMCI)==208 && numel(unique(TMCI.RID))==75, ...
    'Expected 208 MCI observations from 75 participants.');
fprintf('Whole cohort: %d obs / %d participants\n',height(T),numel(unique(T.RID)));
fprintf('MCI subgroup: %d obs / %d participants\n',height(TMCI),numel(unique(TMCI.RID)));

%% Formulas
% Whole cohort: locked primary model (Steps 10/23/30B)
fWhole = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
          'Age + SexModel + Education + ' ...
          'Tau_z:TimeYears + Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
          '(1|RIDcat) + (TimeYears-1|RIDcat)'];

% MCI: mutually adjusted model as in Steps 47/49I (identical fixed effects)
fMCI = ['CDRSB ~ 1 + Age + Education + SexModel + ' ...
        'TimeYears*Tau_z + TimeYears*Centiloid_z + TimeYears*PTau217_z + ' ...
        '(1|RIDcat) + (-1 + TimeYears|RIDcat)'];

%% Fit and test
spec = { ...
    "Whole cohort", T,    fWhole, 'REML'; ...
    "Whole cohort", T,    fWhole, 'ML'; ...
    "MCI",          TMCI, fMCI,   'ML'};

Result = table();
for s = 1:size(spec,1)
    m = fitlme(spec{s,2},spec{s,3},'FitMethod',spec{s,4});
    Result = [Result; localWald(m,spec{s,1},string(spec{s,4}))]; %#ok<AGROW>
end

fprintf('\nRESULTS\n');
disp(Result);

fprintf(['Benchmarks: whole REML diff 0.202 (-0.413 to 0.817), P=.52; ' ...
         'whole ML diff 0.200 (-0.406 to 0.806), P=.52; ' ...
         'MCI ML diff 0.023 (-0.405 to 0.451), P=.92\n']);

%% Save (summary statistics only)
writetable(Result,outFile);
fprintf('\nSaved: %s\n',outFile);
fprintf('============================================================\n');

%% Local function
function row = localWald(m,population,method)
    nm = string(m.CoefficientNames);
    iT = find(nm=="TimeYears:Tau_z"     | nm=="Tau_z:TimeYears");
    iP = find(nm=="TimeYears:PTau217_z" | nm=="PTau217_z:TimeYears");
    assert(isscalar(iT) && isscalar(iP),'Slope terms not uniquely identified.');

    b = m.Coefficients.Estimate;
    V = m.CoefficientCovariance;

    d    = b(iT) - b(iP);
    seD  = sqrt(V(iT,iT) + V(iP,iP) - 2*V(iT,iP));
    rTP  = V(iT,iP) / sqrt(V(iT,iT)*V(iP,iP));
    ci   = d + [-1 1]*tinv(0.975,m.DFE)*seD;

    H = zeros(1,numel(b)); H(iT) = 1; H(iP) = -1;
    [p,F,df1,df2] = coefTest(m,H);

    row = table(population,method,m.NumObservations, ...
        b(iT),sqrt(V(iT,iT)),b(iP),sqrt(V(iP,iP)),rTP, ...
        d,seD,ci(1),ci(2),F,df1,df2,p, ...
        'VariableNames',{'Population','FitMethod','Observations', ...
        'TauTime_Beta','TauTime_SE','PTauTime_Beta','PTauTime_SE', ...
        'SlopeEstimateCorrelation','Difference','Difference_SE', ...
        'Difference_CI_Lower','Difference_CI_Upper','F','DF1','DF2','PValue'});
end
