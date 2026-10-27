% COVARIANCE-STRUCTURE SENSITIVITY ANALYSIS:
% Correlated random-intercept/random-slope model.
% Reported tau x time coefficient: beta = 0.376.
% Estimated intercept-slope correlation approached the
% parameter boundary (+1.00).
% The diagonal covariance structure was retained
% for the primary longitudinal model.

%% Step09_Candidate_Primary_Longitudinal_LME_v3.m
% ADNI Tau PET revised cohort analysis
% Same fitted model as prior Step09 runs.
% Fixes covariance extraction using the MATLAB covarianceParameters(lme)
% function rather than a nonexistent object property.

clear; clc; close all;

outputDir = 'D:\output';
dataFile  = fullfile(outputDir,'Step08_ModelReady_379_Observations.csv');

assert(isfile(dataFile),'Missing Step08 model-ready dataset.');
T = readtable(dataFile,'VariableNamingRule','preserve');

%% Numeric variables
numVars = {'RID','CDRSB','TimeYears','Tau_z','Centiloid_z','PTau217_z','Age','Education'};
for k = 1:numel(numVars)
    if ~isnumeric(T.(numVars{k}))
        T.(numVars{k}) = str2double(string(T.(numVars{k})));
    end
end

%% Factors
T.RIDcat = categorical(T.RID);
T.SexModel = categorical(string(T.Sex),["Male","Female"]);

%% Complete cases
complete = isfinite(T.CDRSB) & isfinite(T.TimeYears) & ...
           isfinite(T.Tau_z) & isfinite(T.Centiloid_z) & ...
           isfinite(T.PTau217_z) & isfinite(T.Age) & ...
           isfinite(T.Education) & ~isundefined(T.SexModel);
TM = T(complete,:);

fprintf('\n====================================================\n');
fprintf(' STEP 9 v3 - CANDIDATE PRIMARY LONGITUDINAL LME\n');
fprintf('====================================================\n');
fprintf('Complete observations: %d\n',height(TM));
fprintf('Complete participants: %d\n',numel(unique(TM.RID)));

formula = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
           'Age + SexModel + Education + Tau_z:TimeYears + ' ...
           'Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
           '(1 + TimeYears|RIDcat)'];

lme = fitlme(TM,formula,'FitMethod','REML');
disp(lme);

%% Fixed effects
FEraw = lme.Coefficients;
if istable(FEraw)
    FE = FEraw;
elseif isa(FEraw,'dataset')
    FE = dataset2table(FEraw);
else
    FE = struct2table(FEraw);
end

names = string(FE.Name);
tauIdx = find(contains(names,'Tau_z') & contains(names,'TimeYears'),1);

tauBeta = FE.Estimate(tauIdx);
tauSE   = FE.SE(tauIdx);
tauP    = FE.pValue(tauIdx);
tauLo   = FE.Lower(tauIdx);
tauHi   = FE.Upper(tauIdx);

fprintf('\nKEY TAU x TIME RESULT\n');
fprintf('Beta = %.6f\n',tauBeta);
fprintf('SE = %.6f\n',tauSE);
fprintf('95%% CI = %.6f to %.6f\n',tauLo,tauHi);
fprintf('P = %.8g\n',tauP);

%% Biomarker x Time summary
terms = ["Tau_z","Centiloid_z","PTau217_z"];
InteractionSummary = table();

for k = 1:numel(terms)
    idx = find(contains(names,terms(k)) & contains(names,'TimeYears'),1);
    if ~isempty(idx)
        row = table(terms(k),FE.Estimate(idx),FE.SE(idx), ...
            FE.Lower(idx),FE.Upper(idx),FE.pValue(idx), ...
            'VariableNames',{'Biomarker','Estimate','SE','CI_Lower','CI_Upper','PValue'});
        InteractionSummary = [InteractionSummary; row]; %#ok<AGROW>
    end
end

%% Correct covariance extraction
[psi,mse,stats] = covarianceParameters(lme);

boundaryFlag = false;
corrEstimate = NaN;

if ~isempty(psi) && ~isempty(psi{1}) && size(psi{1},1) >= 2
    Psi = psi{1};
    sd0 = sqrt(Psi(1,1));
    sd1 = sqrt(Psi(2,2));
    if sd0>0 && sd1>0
        corrEstimate = Psi(1,2)/(sd0*sd1);
        boundaryFlag = abs(corrEstimate) >= 0.999;
    end
else
    Psi = [];
end

fprintf('\nRANDOM-EFFECTS COVARIANCE AUDIT\n');
if ~isempty(Psi)
    fprintf('Random-effects covariance matrix:\n');
    disp(Psi);
end
fprintf('Residual variance (MSE): %.8f\n',mse);
fprintf('Intercept-slope correlation: %.8f\n',corrEstimate);

if boundaryFlag
    fprintf('WARNING: correlation is at/near +/-1 -> boundary/singular fit.\n');
else
    fprintf('No +/-1 boundary correlation detected.\n');
end

%% Model info
ModelInfo = table(height(TM),numel(unique(TM.RID)), ...
    lme.LogLikelihood,lme.ModelCriterion.AIC,lme.ModelCriterion.BIC, ...
    corrEstimate,boundaryFlag, ...
    'VariableNames',{'Observations','Participants','LogLikelihood', ...
                     'AIC','BIC','RandomInterceptSlopeCorrelation','BoundaryFlag'});

%% Save outputs
writetable(FE,fullfile(outputDir,'Step09_PrimaryLME_FixedEffects.csv'));
writetable(InteractionSummary, ...
    fullfile(outputDir,'Step09_PrimaryLME_BiomarkerTime_Interactions.csv'));
writetable(ModelInfo,fullfile(outputDir,'Step09_PrimaryLME_ModelInfo.csv'));

% Save covariance matrix separately
if ~isempty(Psi)
    CovTable = array2table(Psi, ...
        'VariableNames',{'RandomIntercept','RandomTimeSlope'}, ...
        'RowNames',{'RandomIntercept','RandomTimeSlope'});
    writetable(CovTable, ...
        fullfile(outputDir,'Step09_PrimaryLME_RandomEffectsCovariance.csv'), ...
        'WriteRowNames',true);
end

save(fullfile(outputDir,'Step09_PrimaryLME_Model.mat'),'lme','formula','psi','mse','stats');

fprintf('\nOutputs saved successfully to:\n%s\n',outputDir);
fprintf('====================================================\n');
