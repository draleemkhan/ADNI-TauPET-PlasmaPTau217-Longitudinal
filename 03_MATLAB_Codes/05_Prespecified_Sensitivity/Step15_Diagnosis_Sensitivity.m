%% Step15_Diagnosis_Sensitivity.m
% ADNI Tau PET revised cohort analysis
%
% STEP 15 — Diagnosis sensitivity
%
% Fits three REML linear mixed-effects models on the locked
% N=143 / 379-observation dataset, all using:
%   uncorrelated random intercept + random Time slope
%
% Model 1: Primary, no diagnosis
% Model 2: + baseline diagnosis main effect
% Model 3: + baseline diagnosis main effect + diagnosis x Time
%
% Primary fixed effects:
%   CDR-SB ~ Time + Tau + Centiloid + p-tau217 + Age + Sex + Education
%            + Tau*Time + Centiloid*Time + p-tau217*Time
%
% Diagnosis coding:
%   Cognitively normal = reference
%   MCI
%   Dementia
%
% Goal:
% - reproduce Tau x Time beta and P across all three models
% - compare AIC/BIC
% - determine whether Tau x Time remains positive/significant after
%   adjustment for baseline diagnosis and diagnosis-specific progression.
%
% No change to the locked primary model is made by this script.

clear; clc; close all;

outputDir = 'D:\output';

dataFile = fullfile(outputDir,'Step08_ModelReady_379_Observations.csv');
dxFile   = fullfile(outputDir,'Step07_Final143_With_Baseline_Diagnosis.csv');

assert(isfile(dataFile),'Missing Step08 model-ready longitudinal dataset.');
assert(isfile(dxFile),'Missing Step07 diagnosis participant dataset.');

T = readtable(dataFile,'VariableNamingRule','preserve');
D = readtable(dxFile,'VariableNamingRule','preserve');

fprintf('\n====================================================\n');
fprintf(' STEP 15 - DIAGNOSIS SENSITIVITY\n');
fprintf('====================================================\n');

%% Ensure numeric fields
numVars = {'RID','CDRSB','TimeYears','Tau_z','Centiloid_z','PTau217_z','Age','Education'};
for k=1:numel(numVars)
    if ~isnumeric(T.(numVars{k}))
        T.(numVars{k}) = str2double(string(T.(numVars{k})));
    end
end

%% Attach diagnosis if not already present
if ~ismember('Baseline_Diagnosis',T.Properties.VariableNames)
    if ~all(ismember({'RID','Baseline_Diagnosis'},D.Properties.VariableNames))
        error('Diagnosis information unavailable.');
    end
    T = innerjoin(T,D(:,{'RID','Baseline_Diagnosis'}),'Keys','RID');
end

%% Factors
T.RIDcat = categorical(T.RID);
T.SexModel = categorical(string(T.Sex),["Male","Female"]);
T.DiagnosisModel = categorical(string(T.Baseline_Diagnosis), ...
    ["Cognitively normal","MCI","Dementia"]);

%% Complete-case audit
complete = isfinite(T.CDRSB) & isfinite(T.TimeYears) & ...
           isfinite(T.Tau_z) & isfinite(T.Centiloid_z) & ...
           isfinite(T.PTau217_z) & isfinite(T.Age) & ...
           isfinite(T.Education) & ~isundefined(T.SexModel) & ...
           ~isundefined(T.DiagnosisModel);

M = T(complete,:);

fprintf('Complete observations: %d [expected 379]\n',height(M));
fprintf('Complete participants: %d [expected 143]\n',numel(unique(M.RID)));

fprintf('\nDIAGNOSIS COUNTS\n');
cats = categories(M.DiagnosisModel);
for k=1:numel(cats)
    fprintf('%s: %d participants\n',cats{k}, ...
        numel(unique(M.RID(M.DiagnosisModel==cats{k}))));
end

%% Model formulas
f1 = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
      'Age + SexModel + Education + Tau_z:TimeYears + ' ...
      'Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
      '(1 + TimeYears|RIDcat)'];

f2 = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
      'Age + SexModel + Education + DiagnosisModel + ' ...
      'Tau_z:TimeYears + Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
      '(1 + TimeYears|RIDcat)'];

f3 = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
      'Age + SexModel + Education + DiagnosisModel + ' ...
      'Tau_z:TimeYears + Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
      'DiagnosisModel:TimeYears + (1 + TimeYears|RIDcat)'];

%% Fit all three with diagonal covariance
fprintf('\nFitting Model 1: Primary, no diagnosis...\n');
m1 = fitlme(M,f1,'FitMethod','REML','CovariancePattern','Diagonal');

fprintf('Fitting Model 2: + diagnosis main effect...\n');
m2 = fitlme(M,f2,'FitMethod','REML','CovariancePattern','Diagonal');

fprintf('Fitting Model 3: + diagnosis and diagnosis x Time...\n');
m3 = fitlme(M,f3,'FitMethod','REML','CovariancePattern','Diagonal');

%% Extract summaries
models = {m1,m2,m3};
labels = ["Primary: no diagnosis"; ...
          "+ diagnosis main effect"; ...
          "+ diagnosis and diagnosis x Time"];

Result = table();

for k=1:3
    FEraw = models{k}.Coefficients;
    if istable(FEraw)
        FE = FEraw;
    elseif isa(FEraw,'dataset')
        FE = dataset2table(FEraw);
    else
        try
            FE = dataset2table(FEraw);
        catch
            FE = struct2table(FEraw);
        end
    end

    names = string(FE.Name);
    ix = find(names=="TimeYears:Tau_z" | names=="Tau_z:TimeYears",1);
    if isempty(ix)
        error('Tau x Time term not found in model %d.',k);
    end

    row = table(labels(k), ...
        FE.Estimate(ix),FE.SE(ix),FE.Lower(ix),FE.Upper(ix),FE.pValue(ix), ...
        models{k}.ModelCriterion.AIC, ...
        models{k}.ModelCriterion.BIC, ...
        models{k}.LogLikelihood, ...
        'VariableNames',{'Model','TauTime_Beta','SE','CI_Lower','CI_Upper', ...
                         'PValue','AIC','BIC','LogLikelihood'});
    Result = [Result; row]; %#ok<AGROW>
end

fprintf('\n====================================================\n');
fprintf(' DIAGNOSIS SENSITIVITY SUMMARY\n');
fprintf('====================================================\n');
disp(Result);

%% Extract diagnosis terms from Models 2 and 3
FE2 = convertFE(m2.Coefficients);
FE3 = convertFE(m3.Coefficients);

DxTerms2 = FE2(contains(string(FE2.Name),'DiagnosisModel'),:);
DxTerms3 = FE3(contains(string(FE3.Name),'DiagnosisModel'),:);

fprintf('\nMODEL 2 DIAGNOSIS TERMS\n');
disp(DxTerms2);

fprintf('\nMODEL 3 DIAGNOSIS / DIAGNOSIS x TIME TERMS\n');
disp(DxTerms3);

%% Compare AIC/BIC changes relative to primary
Result.DeltaAIC_vsPrimary = Result.AIC - Result.AIC(1);
Result.DeltaBIC_vsPrimary = Result.BIC - Result.BIC(1);

fprintf('\nWITH DELTA INFORMATION CRITERIA\n');
disp(Result);

%% Direction and significance audit
Result.PositiveTauTime = Result.TauTime_Beta > 0;
Result.Significant_0p05 = Result.PValue < 0.05;

fprintf('\nTau x Time positive in all models: %s\n', ...
    string(all(Result.PositiveTauTime)));
fprintf('Tau x Time P<0.05 in all models: %s\n', ...
    string(all(Result.Significant_0p05)));

%% Save outputs
writetable(Result, ...
    fullfile(outputDir,'Step15_Diagnosis_Sensitivity_ModelComparison.csv'));

writetable(FE2, ...
    fullfile(outputDir,'Step15_Model2_DiagnosisMain_FixedEffects.csv'));

writetable(FE3, ...
    fullfile(outputDir,'Step15_Model3_DiagnosisTime_FixedEffects.csv'));

save(fullfile(outputDir,'Step15_Diagnosis_Sensitivity_Models.mat'), ...
    'm1','m2','m3','f1','f2','f3');

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');

%% Local helper
function FE = convertFE(FEraw)
    if istable(FEraw)
        FE = FEraw;
    elseif isa(FEraw,'dataset')
        FE = dataset2table(FEraw);
    else
        try
            FE = dataset2table(FEraw);
        catch
            FE = struct2table(FEraw);
        end
    end
end
