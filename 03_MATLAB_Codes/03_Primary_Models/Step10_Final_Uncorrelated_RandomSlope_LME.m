%% Step10_Final_Uncorrelated_RandomSlope_LME.m
% ADNI Tau PET revised cohort analysis
%
% LOCKED FINAL PRIMARY RANDOM-EFFECTS STRUCTURE:
%   Uncorrelated random intercept + random Time slope
%
% Fixed effects:
%   CDR-SB ~ Time + Tau + Centiloid + p-tau217 + Age + Sex + Education
%            + Tau*Time + Centiloid*Time + p-tau217*Time
%
% Random effects:
%   (1 + TimeYears | RIDcat), covariance constrained to DIAGONAL
%   so intercept-slope covariance = 0.
%
% Fit method:
%   REML
%
% Key parameter:
%   Tau_z:TimeYears

clear; clc; close all;

%% Locked path
outputDir = 'D:\output';
dataFile  = fullfile(outputDir,'Step08_ModelReady_379_Observations.csv');

assert(isfile(dataFile),'Missing Step08 model-ready dataset.');

T = readtable(dataFile,'VariableNamingRule','preserve');

fprintf('\n====================================================\n');
fprintf(' STEP 10 - FINAL UNCORRELATED RANDOM-SLOPE LME\n');
fprintf('====================================================\n');

%% Ensure numeric fields
numVars = {'RID','CDRSB','TimeYears','Tau_z','Centiloid_z','PTau217_z', ...
           'Age','Education'};

for k = 1:numel(numVars)
    v = numVars{k};
    if ~isnumeric(T.(v))
        T.(v) = str2double(string(T.(v)));
    end
end

%% Factors
T.RIDcat = categorical(T.RID);
T.SexModel = categorical(string(T.Sex),["Male","Female"]);

%% Complete-case audit
complete = isfinite(T.CDRSB) & isfinite(T.TimeYears) & ...
           isfinite(T.Tau_z) & isfinite(T.Centiloid_z) & ...
           isfinite(T.PTau217_z) & isfinite(T.Age) & ...
           isfinite(T.Education) & ~isundefined(T.SexModel);

TM = T(complete,:);

fprintf('Complete observations: %d [expected 379]\n',height(TM));
fprintf('Complete participants: %d [expected 143]\n',numel(unique(TM.RID)));

%% Locked formula
formula = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
           'Age + SexModel + Education + Tau_z:TimeYears + ' ...
           'Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
           '(1 + TimeYears|RIDcat)'];

fprintf('\nLOCKED FORMULA\n%s\n',formula);
fprintf('Covariance pattern: Diagonal (intercept-slope covariance fixed at 0)\n');

%% Fit final primary model
lme = fitlme(TM,formula, ...
    'FitMethod','REML', ...
    'CovariancePattern','Diagonal');

fprintf('\nFINAL PRIMARY MODEL FIT COMPLETE\n');
disp(lme);

%% Fixed effects
FEraw = lme.Coefficients;

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

fprintf('\nFIXED EFFECTS\n');
disp(FE);

%% Key Tau x Time coefficient
names = string(FE.Name);
tauIdx = find(contains(names,'Tau_z') & contains(names,'TimeYears'),1);

if isempty(tauIdx)
    error('Tau x Time coefficient not found.');
end

tauBeta = FE.Estimate(tauIdx);
tauSE   = FE.SE(tauIdx);
tauLo   = FE.Lower(tauIdx);
tauHi   = FE.Upper(tauIdx);
tauP    = FE.pValue(tauIdx);
tauT    = FE.tStat(tauIdx);
tauDF   = FE.DF(tauIdx);

fprintf('\n====================================================\n');
fprintf(' FINAL PRIMARY RESULT: TAU x TIME\n');
fprintf('====================================================\n');
fprintf('Beta = %.6f CDR-SB points/year per 1-SD higher tau\n',tauBeta);
fprintf('SE = %.6f\n',tauSE);
fprintf('95%% CI = %.6f to %.6f\n',tauLo,tauHi);
fprintf('t = %.4f | df = %.2f | P = %.8g\n',tauT,tauDF,tauP);

%% Biomarker x Time interactions
terms = ["Tau_z","Centiloid_z","PTau217_z"];
InteractionSummary = table();

for k = 1:numel(terms)
    idx = find(contains(names,terms(k)) & contains(names,'TimeYears'),1);
    if ~isempty(idx)
        row = table(terms(k),FE.Estimate(idx),FE.SE(idx), ...
            FE.Lower(idx),FE.Upper(idx),FE.pValue(idx), ...
            'VariableNames',{'Biomarker','Estimate','SE', ...
                             'CI_Lower','CI_Upper','PValue'});
        InteractionSummary = [InteractionSummary; row]; %#ok<AGROW>
    end
end

fprintf('\nBIOMARKER x TIME INTERACTIONS\n');
disp(InteractionSummary);

%% Covariance parameters
[psi,mse,stats] = covarianceParameters(lme);

fprintf('\nRANDOM-EFFECTS COVARIANCE AUDIT\n');

if ~isempty(psi) && ~isempty(psi{1})
    Psi = psi{1};
    disp(Psi);

    if size(Psi,1)>=2
        cov01 = Psi(1,2);
        fprintf('Estimated intercept-slope covariance: %.10f\n',cov01);

        if abs(cov01) < 1e-10
            fprintf('UNCORRELATED RANDOM-EFFECTS CONSTRAINT: PASS\n');
        else
            fprintf('UNCORRELATED RANDOM-EFFECTS CONSTRAINT: REVIEW\n');
        end

        fprintf('Random intercept SD: %.6f\n',sqrt(Psi(1,1)));
        fprintf('Random Time slope SD: %.6f\n',sqrt(Psi(2,2)));
    end
else
    Psi = [];
end

fprintf('Residual SD: %.6f\n',sqrt(mse));

%% Model information
ModelInfo = table( ...
    height(TM), ...
    numel(unique(TM.RID)), ...
    lme.LogLikelihood, ...
    lme.ModelCriterion.AIC, ...
    lme.ModelCriterion.BIC, ...
    tauBeta,tauSE,tauLo,tauHi,tauP, ...
    'VariableNames',{ ...
    'Observations','Participants','LogLikelihood','AIC','BIC', ...
    'TauTime_Beta','TauTime_SE','TauTime_CI_Lower', ...
    'TauTime_CI_Upper','TauTime_PValue'});

%% Save outputs
writetable(FE, ...
    fullfile(outputDir,'Step10_FinalPrimary_UncorrelatedLME_FixedEffects.csv'));

writetable(InteractionSummary, ...
    fullfile(outputDir,'Step10_FinalPrimary_BiomarkerTime_Interactions.csv'));

writetable(ModelInfo, ...
    fullfile(outputDir,'Step10_FinalPrimary_ModelInfo.csv'));

if ~isempty(Psi)
    CovTable = array2table(Psi, ...
        'VariableNames',{'RandomIntercept','RandomTimeSlope'}, ...
        'RowNames',{'RandomIntercept','RandomTimeSlope'});

    writetable(CovTable, ...
        fullfile(outputDir,'Step10_FinalPrimary_RandomEffectsCovariance.csv'), ...
        'WriteRowNames',true);
end

save(fullfile(outputDir,'Step10_FinalPrimary_UncorrelatedLME_Model.mat'), ...
    'lme','formula','psi','mse','stats');

%% Text report
txtFile = fullfile(outputDir,'Step10_FinalPrimary_UncorrelatedLME_Report.txt');
fid = fopen(txtFile,'w');

fprintf(fid,'FINAL PRIMARY ADNI TAU PET LONGITUDINAL MODEL\n');
fprintf(fid,'============================================\n\n');
fprintf(fid,'Participants: %d\n',numel(unique(TM.RID)));
fprintf(fid,'Observations: %d\n',height(TM));
fprintf(fid,'Fit method: REML\n');
fprintf(fid,'Random effects: uncorrelated random intercept + random Time slope\n');
fprintf(fid,'Covariance pattern: Diagonal\n\n');

fprintf(fid,'TAU x TIME PRIMARY RESULT\n');
fprintf(fid,'Beta = %.8f\n',tauBeta);
fprintf(fid,'SE = %.8f\n',tauSE);
fprintf(fid,'95%% CI = %.8f to %.8f\n',tauLo,tauHi);
fprintf(fid,'t = %.8f\n',tauT);
fprintf(fid,'df = %.3f\n',tauDF);
fprintf(fid,'P = %.10g\n\n',tauP);

fprintf(fid,'AIC = %.6f\n',lme.ModelCriterion.AIC);
fprintf(fid,'BIC = %.6f\n',lme.ModelCriterion.BIC);
fprintf(fid,'LogLikelihood = %.6f\n',lme.LogLikelihood);

fclose(fid);

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');
