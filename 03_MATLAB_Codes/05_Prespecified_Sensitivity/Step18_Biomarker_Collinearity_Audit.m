%% Step18_Biomarker_Collinearity_Audit.m
% ADNI Tau PET revised cohort analysis
%
% STEP 18 — Biomarker collinearity audit
%
% PURPOSE
% Quantify correlation among baseline biomarkers and multivariable
% collinearity within the locked primary longitudinal model.
%
% Correlations are computed at the PARTICIPANT level (N=143) using the
% baseline values:
%   Tau PET meta-temporal SUVR
%   Centiloid
%   Plasma p-tau217
%
% VIFs are computed from the FIXED-EFFECT design matrix of the locked
% primary longitudinal model across all 379 repeated observations:
%
%   Time + Tau + Centiloid + p-tau217 + Age + Sex + Education
%   + Tau*Time + Centiloid*Time + p-tau217*Time
%
% Historical approximate benchmarks:
%   Tau PET vs p-tau217 Pearson r ~0.81
%   Tau PET vs Centiloid Pearson r ~0.51
%   Centiloid vs p-tau217 Pearson r ~0.51
%   Tau PET VIF ~5.4
%   p-tau217 VIF ~5.5
%   Centiloid VIF ~2.7
%   Tau x Time VIF ~4.6
%   p-tau217 x Time VIF ~4.8
%
% This is a diagnostic audit only and does not modify the primary model.

clear; clc; close all;

outputDir = 'D:\output';

participantFile = fullfile(outputDir,'Step08_Final143_With_Demographics.csv');
modelFile       = fullfile(outputDir,'Step08_ModelReady_379_Observations.csv');

assert(isfile(participantFile),'Missing Step08 participant-level file.');
assert(isfile(modelFile),'Missing Step08 model-ready longitudinal file.');

P = readtable(participantFile,'VariableNamingRule','preserve');
T = readtable(modelFile,'VariableNamingRule','preserve');

fprintf('\n====================================================\n');
fprintf(' STEP 18 - BIOMARKER COLLINEARITY AUDIT\n');
fprintf('====================================================\n');

%% Participant-level biomarker variables
tau = double(P.META_TEMPORAL_SUVR);
cent = double(P.Centiloid);
ptau = double(P.PTau217_Fujirebio);

ok = isfinite(tau) & isfinite(cent) & isfinite(ptau);

fprintf('Participants with complete three-biomarker data: %d / %d\n',sum(ok),height(P));

%% Pearson correlations
[rTauPtau,pTauPtau] = corr(tau(ok),ptau(ok),'Type','Pearson');
[rTauCent,pTauCent] = corr(tau(ok),cent(ok),'Type','Pearson');
[rCentPtau,pCentPtau] = corr(cent(ok),ptau(ok),'Type','Pearson');

fprintf('\nPARTICIPANT-LEVEL PEARSON CORRELATIONS\n');
fprintf('Tau PET vs p-tau217:   r = %.4f, P = %.6g\n',rTauPtau,pTauPtau);
fprintf('Tau PET vs Centiloid:  r = %.4f, P = %.6g\n',rTauCent,pTauCent);
fprintf('Centiloid vs p-tau217: r = %.4f, P = %.6g\n',rCentPtau,pCentPtau);

CorrelationTable = table( ...
    ["Tau PET vs p-tau217";"Tau PET vs Centiloid";"Centiloid vs p-tau217"], ...
    [rTauPtau;rTauCent;rCentPtau], ...
    [pTauPtau;pTauCent;pCentPtau], ...
    'VariableNames',{'Comparison','PearsonR','PValue'});

%% Prepare fixed-effect design matrix
numVars = {'TimeYears','Tau_z','Centiloid_z','PTau217_z','Age','Education'};
for k=1:numel(numVars)
    if ~isnumeric(T.(numVars{k}))
        T.(numVars{k}) = str2double(string(T.(numVars{k})));
    end
end

sex = string(T.Sex);
Female = double(sex=="Female");

complete = isfinite(T.TimeYears) & isfinite(T.Tau_z) & ...
           isfinite(T.Centiloid_z) & isfinite(T.PTau217_z) & ...
           isfinite(T.Age) & isfinite(T.Education) & ...
           ismember(sex,["Male","Female"]);

M = T(complete,:);
Female = Female(complete);

fprintf('\nFixed-effect design rows for VIF: %d [expected 379]\n',height(M));

%% Fixed-effect predictors, excluding intercept
X = [ ...
    M.TimeYears, ...
    M.Tau_z, ...
    M.Centiloid_z, ...
    M.PTau217_z, ...
    M.Age, ...
    Female, ...
    M.Education, ...
    M.TimeYears .* M.Tau_z, ...
    M.TimeYears .* M.Centiloid_z, ...
    M.TimeYears .* M.PTau217_z];

names = [ ...
    "TimeYears"
    "Tau_z"
    "Centiloid_z"
    "PTau217_z"
    "Age"
    "Female"
    "Education"
    "Tau_z x Time"
    "Centiloid_z x Time"
    "PTau217_z x Time"];

%% Calculate VIF for each predictor
p = size(X,2);
VIF = nan(p,1);
Tolerance = nan(p,1);
Rsq = nan(p,1);

for j=1:p
    yj = X(:,j);
    Xother = X(:,setdiff(1:p,j));

    % Add intercept for auxiliary regression.
    Xaux = [ones(size(Xother,1),1), Xother];
    b = Xaux \ yj;
    yhat = Xaux*b;

    SSE = sum((yj-yhat).^2);
    SST = sum((yj-mean(yj)).^2);

    if SST <= eps
        Rsq(j)=NaN;
        Tolerance(j)=NaN;
        VIF(j)=Inf;
    else
        Rsq(j)=1-SSE/SST;
        Tolerance(j)=1-Rsq(j);
        VIF(j)=1/Tolerance(j);
    end
end

VIFTable = table(names,VIF,Tolerance,Rsq, ...
    'VariableNames',{'Predictor','VIF','Tolerance','Auxiliary_R2'});

fprintf('\nFIXED-EFFECT VIF AUDIT\n');
disp(VIFTable);

%% Key benchmark extract
keyNames = ["Tau_z","PTau217_z","Centiloid_z","Tau_z x Time","PTau217_z x Time","Centiloid_z x Time"];
KeyVIF = VIFTable(ismember(VIFTable.Predictor,keyNames),:);

fprintf('\nKEY BIOMARKER VIFS\n');
disp(KeyVIF);

%% Condition diagnostics
% Standardize continuous predictor columns before condition number so units
% do not dominate. Female remains binary and is standardized likewise.
Xs = zscore(X);
[~,S,~] = svd(Xs,'econ');
singVals = diag(S);
conditionNumber = max(singVals)/min(singVals);

fprintf('\nDESIGN-MATRIX CONDITION NUMBER (standardized predictors): %.4f\n',conditionNumber);

%% Save outputs
writetable(CorrelationTable, ...
    fullfile(outputDir,'Step18_Biomarker_Correlations.csv'));

writetable(VIFTable, ...
    fullfile(outputDir,'Step18_PrimaryModel_VIF_AllPredictors.csv'));

writetable(KeyVIF, ...
    fullfile(outputDir,'Step18_PrimaryModel_VIF_KeyBiomarkers.csv'));

Summary = table( ...
    sum(ok),rTauPtau,rTauCent,rCentPtau,conditionNumber, ...
    'VariableNames',{'Participants','Tau_vs_PTau217_r','Tau_vs_Centiloid_r', ...
                     'Centiloid_vs_PTau217_r','StandardizedDesignConditionNumber'});

writetable(Summary, ...
    fullfile(outputDir,'Step18_Biomarker_Collinearity_Summary.csv'));

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');
