%% Step20_Incremental_Biomarker_Models_ML.m
% ADNI Tau PET revised cohort analysis
%
% STEP 20 — Incremental biomarker models using MAXIMUM LIKELIHOOD
%
% PURPOSE
% Formally evaluate incremental longitudinal model fit across nested
% fixed-effect specifications while keeping the same locked random-effects
% structure:
%   uncorrelated random intercept + random Time slope
%
% Model 1:
%   Centiloid + demographics
%
% Model 2:
%   Model 1 + plasma p-tau217
%
% Model 3:
%   Model 2 + tau PET
%
% Each biomarker enters as BOTH a main effect and a Time interaction.
%
% ML is used because AIC/BIC and likelihood-ratio comparisons across
% different fixed-effect specifications should be based on ML fits.
%
% Outputs:
%   - exact AIC/BIC/logLik
%   - key biomarker x Time coefficient for each model
%   - nested likelihood-ratio tests:
%       Model 2 vs Model 1
%       Model 3 vs Model 2
%
% This step complements Step19 marginal/conditional R2.

clear; clc; close all;

outputDir = 'D:\output';
dataFile = fullfile(outputDir,'Step08_ModelReady_379_Observations.csv');

assert(isfile(dataFile),'Missing Step08 model-ready dataset.');

T = readtable(dataFile,'VariableNamingRule','preserve');

fprintf('\n====================================================\n');
fprintf(' STEP 20 - INCREMENTAL BIOMARKER MODELS (ML)\n');
fprintf('====================================================\n');

%% Ensure numeric fields
numVars = {'RID','CDRSB','TimeYears','Tau_z','Centiloid_z','PTau217_z','Age','Education'};
for k=1:numel(numVars)
    if ~isnumeric(T.(numVars{k}))
        T.(numVars{k}) = str2double(string(T.(numVars{k})));
    end
end

T.RIDcat = categorical(T.RID);
T.SexModel = categorical(string(T.Sex),["Male","Female"]);

complete = isfinite(T.CDRSB) & isfinite(T.TimeYears) & ...
           isfinite(T.Tau_z) & isfinite(T.Centiloid_z) & ...
           isfinite(T.PTau217_z) & isfinite(T.Age) & ...
           isfinite(T.Education) & ~isundefined(T.SexModel);

M = T(complete,:);

fprintf('Complete observations: %d [expected 379]\n',height(M));
fprintf('Complete participants: %d [expected 143]\n',numel(unique(M.RID)));

%% Nested fixed-effect formulas
f1 = ['CDRSB ~ TimeYears + Centiloid_z + Age + SexModel + Education + ' ...
      'Centiloid_z:TimeYears + (1 + TimeYears|RIDcat)'];

f2 = ['CDRSB ~ TimeYears + Centiloid_z + PTau217_z + Age + SexModel + Education + ' ...
      'Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
      '(1 + TimeYears|RIDcat)'];

f3 = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + Age + SexModel + Education + ' ...
      'Tau_z:TimeYears + Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
      '(1 + TimeYears|RIDcat)'];

labels = ["Centiloid + demographics"; ...
          "+ p-tau217"; ...
          "+ tau PET"];

%% Fit with ML
fprintf('\nFitting Model 1 with ML...\n');
m1 = fitlme(M,f1,'FitMethod','ML','CovariancePattern','Diagonal');

fprintf('Fitting Model 2 with ML...\n');
m2 = fitlme(M,f2,'FitMethod','ML','CovariancePattern','Diagonal');

fprintf('Fitting Model 3 with ML...\n');
m3 = fitlme(M,f3,'FitMethod','ML','CovariancePattern','Diagonal');

models = {m1,m2,m3};

%% Extract key interaction from each model
Result = table();

for k=1:3
    FEraw = models{k}.Coefficients;
    if istable(FEraw)
        FE = FEraw;
    else
        FE = dataset2table(FEraw);
    end

    names = string(FE.Name);

    if k==1
        ix = find(names=="TimeYears:Centiloid_z" | names=="Centiloid_z:TimeYears",1);
        keyName = "Centiloid x Time";
    elseif k==2
        ix = find(names=="TimeYears:PTau217_z" | names=="PTau217_z:TimeYears",1);
        keyName = "p-tau217 x Time";
    else
        ix = find(names=="TimeYears:Tau_z" | names=="Tau_z:TimeYears",1);
        keyName = "Tau PET x Time";
    end

    if isempty(ix)
        error('Key biomarker x Time term not found in model %d.',k);
    end

    row = table(labels(k),keyName, ...
        FE.Estimate(ix),FE.SE(ix),FE.Lower(ix),FE.Upper(ix),FE.pValue(ix), ...
        models{k}.ModelCriterion.AIC, ...
        models{k}.ModelCriterion.BIC, ...
        models{k}.LogLikelihood, ...
        'VariableNames',{'Model','KeyInteraction','Beta','SE','CI_Lower', ...
                         'CI_Upper','PValue','AIC','BIC','LogLikelihood'});

    Result = [Result;row]; %#ok<AGROW>
end

%% Nested likelihood-ratio tests
% MATLAB compare returns model-comparison statistics for nested LME models.
cmp21 = compare(m1,m2,'CheckNesting',true);
cmp32 = compare(m2,m3,'CheckNesting',true);

fprintf('\n====================================================\n');
fprintf(' INCREMENTAL MODEL SUMMARY\n');
fprintf('====================================================\n');
disp(Result);

fprintf('\nMODEL 2 vs MODEL 1 (add p-tau217)\n');
disp(cmp21);

fprintf('\nMODEL 3 vs MODEL 2 (add tau PET)\n');
disp(cmp32);

%% Parse comparison p-values robustly
p21 = NaN; p32 = NaN;
try
    if istable(cmp21)
        if ismember('pValue',cmp21.Properties.VariableNames)
            p21 = cmp21.pValue(end);
        end
    elseif isa(cmp21,'dataset')
        c21 = dataset2table(cmp21);
        if ismember('pValue',c21.Properties.VariableNames)
            p21 = c21.pValue(end);
        end
    end
catch
end

try
    if istable(cmp32)
        if ismember('pValue',cmp32.Properties.VariableNames)
            p32 = cmp32.pValue(end);
        end
    elseif isa(cmp32,'dataset')
        c32 = dataset2table(cmp32);
        if ismember('pValue',c32.Properties.VariableNames)
            p32 = c32.pValue(end);
        end
    end
catch
end

%% Delta information criteria
Result.DeltaAIC_vsPrevious = [NaN; diff(Result.AIC)];
Result.DeltaBIC_vsPrevious = [NaN; diff(Result.BIC)];

fprintf('\nWITH DELTA INFORMATION CRITERIA\n');
disp(Result);

fprintf('\nINCREMENTAL FIT TESTS\n');
fprintf('Add p-tau217: LRT P = %.8g\n',p21);
fprintf('Add tau PET:  LRT P = %.8g\n',p32);

%% Save outputs
writetable(Result, ...
    fullfile(outputDir,'Step20_Incremental_Biomarker_Models_ML.csv'));

% Save comparison objects as tables where possible
if isa(cmp21,'dataset')
    c21 = dataset2table(cmp21);
else
    c21 = cmp21;
end
if isa(cmp32,'dataset')
    c32 = dataset2table(cmp32);
else
    c32 = cmp32;
end

if istable(c21)
    writetable(c21,fullfile(outputDir,'Step20_LRT_Add_PTau217.csv'));
end
if istable(c32)
    writetable(c32,fullfile(outputDir,'Step20_LRT_Add_TauPET.csv'));
end

Summary = table( ...
    Result.AIC(1),Result.BIC(1), ...
    Result.AIC(2),Result.BIC(2), ...
    Result.AIC(3),Result.BIC(3), ...
    p21,p32, ...
    'VariableNames',{'Model1_AIC','Model1_BIC','Model2_AIC','Model2_BIC', ...
    'Model3_AIC','Model3_BIC','LRT_Add_PTau217_P','LRT_Add_TauPET_P'});

writetable(Summary, ...
    fullfile(outputDir,'Step20_Incremental_Biomarker_Models_ML_Summary.csv'));

save(fullfile(outputDir,'Step20_Incremental_Biomarker_Models_ML.mat'), ...
    'm1','m2','m3','f1','f2','f3','Result','cmp21','cmp32');

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');
