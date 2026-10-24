%% Step13_Residual_and_Influence_Diagnostics.m
% ADNI Tau PET revised cohort analysis
%
% PURPOSE
% Residual and influence diagnostics for the LOCKED primary LMM:
%   uncorrelated random intercept + random Time slope
%
% Benchmarks to reproduce:
% - residuals centered with limited ordinary skewness
% - heavy tails and heteroscedasticity vs fitted CDR-SB
% - only ~2/379 observations with |standardized residual| > 3
% - removing strongest residual subjects individually or jointly does not
%   eliminate the Tau x Time association
% - historical joint-removal benchmark ~ beta 0.38, P ~ .014
%
% IMPORTANT
% No observation is excluded from the primary analysis based on diagnostics.

clear; clc; close all;

outputDir = 'D:\output';
dataFile  = fullfile(outputDir,'Step08_ModelReady_379_Observations.csv');

assert(isfile(dataFile),'Missing Step08 model-ready dataset.');

T = readtable(dataFile,'VariableNamingRule','preserve');

fprintf('\n====================================================\n');
fprintf(' STEP 13 - RESIDUAL AND INFLUENCE DIAGNOSTICS\n');
fprintf('====================================================\n');

%% Ensure variables
numVars = {'RID','CDRSB','TimeYears','Tau_z','Centiloid_z','PTau217_z','Age','Education'};
for k=1:numel(numVars)
    if ~isnumeric(T.(numVars{k}))
        T.(numVars{k})=str2double(string(T.(numVars{k})));
    end
end

T.RIDcat = categorical(T.RID);
T.SexModel = categorical(string(T.Sex),["Male","Female"]);

complete = isfinite(T.CDRSB) & isfinite(T.TimeYears) & ...
           isfinite(T.Tau_z) & isfinite(T.Centiloid_z) & ...
           isfinite(T.PTau217_z) & isfinite(T.Age) & ...
           isfinite(T.Education) & ~isundefined(T.SexModel);

M = T(complete,:);

fprintf('Complete observations: %d\n',height(M));
fprintf('Complete participants: %d\n',numel(unique(M.RID)));

%% Locked primary model
formula = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
           'Age + SexModel + Education + Tau_z:TimeYears + ' ...
           'Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
           '(1 + TimeYears|RIDcat)'];

lme = fitlme(M,formula,'FitMethod','REML','CovariancePattern','Diagonal');

%% Fitted values and residuals
fitted = fitted(lme);
residRaw = residuals(lme,'ResidualType','Raw');

% Pearson-like standardized residual using residual SD
[~,mse] = covarianceParameters(lme);
resSD = sqrt(mse);
stdResid = residRaw ./ resSD;

%% Residual summaries
rMean = mean(residRaw,'omitnan');
rMedian = median(residRaw,'omitnan');
rSD = std(residRaw,'omitnan');
rSkew = skewness(residRaw,0);
rKurt = kurtosis(residRaw,0);

nAbs3 = sum(abs(stdResid)>3);

fprintf('\nRESIDUAL SUMMARY\n');
fprintf('Mean raw residual: %.6f\n',rMean);
fprintf('Median raw residual: %.6f\n',rMedian);
fprintf('Residual SD: %.6f\n',rSD);
fprintf('Skewness: %.4f\n',rSkew);
fprintf('Kurtosis: %.4f\n',rKurt);
fprintf('|standardized residual| > 3: %d / %d\n',nAbs3,height(M));

%% Heteroscedasticity audit
% Correlation between absolute residual magnitude and fitted value
[rho,pRho] = corr(abs(stdResid),fitted,'Type','Spearman','Rows','complete');

% Also compare residual spread across fitted-value quartiles
qEdges = quantile(fitted,[0 .25 .5 .75 1]);
qGroup = discretize(fitted,qEdges,'IncludedEdge','right');
spread = nan(4,1);
for q=1:4
    spread(q)=std(residRaw(qGroup==q),'omitnan');
end

fprintf('\nHETEROSCEDASTICITY AUDIT\n');
fprintf('Spearman rho(|std residual|, fitted) = %.4f, P = %.6g\n',rho,pRho);
fprintf('Residual SD by fitted quartile: %.4f | %.4f | %.4f | %.4f\n',spread);

%% Observation-level extreme residual audit
ObsAudit = table(M.RID,M.TimeYears,M.CDRSB,fitted,residRaw,stdResid, ...
    'VariableNames',{'RID','TimeYears','Observed_CDRSB','Fitted_CDRSB', ...
                     'RawResidual','StandardizedResidual'});

ObsAudit.AbsStdResidual = abs(ObsAudit.StandardizedResidual);
ObsAudit = sortrows(ObsAudit,'AbsStdResidual','descend');

fprintf('\nTOP 10 OBSERVATION-LEVEL RESIDUALS\n');
disp(ObsAudit(1:min(10,height(ObsAudit)),:));

%% Participant-level residual influence ranking
uRID = unique(M.RID);
nSub = numel(uRID);

SubAudit = table();
SubAudit.RID = uRID;
SubAudit.Nobs = zeros(nSub,1);
SubAudit.MaxAbsStdResidual = nan(nSub,1);
SubAudit.MeanAbsStdResidual = nan(nSub,1);
SubAudit.SumSqStdResidual = nan(nSub,1);

for i=1:nSub
    g = M.RID==uRID(i);
    rr = stdResid(g);
    SubAudit.Nobs(i)=sum(g);
    SubAudit.MaxAbsStdResidual(i)=max(abs(rr));
    SubAudit.MeanAbsStdResidual(i)=mean(abs(rr));
    SubAudit.SumSqStdResidual(i)=sum(rr.^2);
end

SubAudit = sortrows(SubAudit,{'MaxAbsStdResidual','SumSqStdResidual'}, ...
    {'descend','descend'});

fprintf('\nTOP 10 PARTICIPANTS BY MAX |STANDARDIZED RESIDUAL|\n');
disp(SubAudit(1:min(10,height(SubAudit)),:));

%% Primary Tau x Time result
FE = dataset2table(lme.Coefficients);
names = string(FE.Name);
tauIdx = find(names=="TimeYears:Tau_z" | names=="Tau_z:TimeYears",1);

primaryBeta = FE.Estimate(tauIdx);
primaryP = FE.pValue(tauIdx);
primaryLo = FE.Lower(tauIdx);
primaryHi = FE.Upper(tauIdx);

fprintf('\nPRIMARY MODEL TAU x TIME\n');
fprintf('Beta = %.6f | 95%% CI %.6f to %.6f | P = %.8g\n', ...
    primaryBeta,primaryLo,primaryHi,primaryP);

%% Leave-out strongest residual subjects
topRID = SubAudit.RID(1:min(2,height(SubAudit)));

InfluenceResults = table();

% Remove each top subject individually
for j=1:numel(topRID)
    keep = M.RID~=topRID(j);
    X = M(keep,:);
    mdl = fitlme(X,formula,'FitMethod','REML','CovariancePattern','Diagonal');
    FEx = dataset2table(mdl.Coefficients);
    nn = string(FEx.Name);
    ix = find(nn=="TimeYears:Tau_z" | nn=="Tau_z:TimeYears",1);

    row = table("Remove RID "+string(topRID(j)), ...
        topRID(j),height(X),numel(unique(X.RID)), ...
        FEx.Estimate(ix),FEx.SE(ix),FEx.Lower(ix),FEx.Upper(ix),FEx.pValue(ix), ...
        'VariableNames',{'Analysis','RemovedRID','Observations','Participants', ...
                         'TauTime_Beta','SE','CI_Lower','CI_Upper','PValue'});
    InfluenceResults = [InfluenceResults; row]; %#ok<AGROW>
end

% Remove both top subjects together
keep = ~ismember(M.RID,topRID);
X = M(keep,:);
mdl = fitlme(X,formula,'FitMethod','REML','CovariancePattern','Diagonal');
FEx = dataset2table(mdl.Coefficients);
nn = string(FEx.Name);
ix = find(nn=="TimeYears:Tau_z" | nn=="Tau_z:TimeYears",1);

row = table("Remove top 2 residual subjects", ...
    NaN,height(X),numel(unique(X.RID)), ...
    FEx.Estimate(ix),FEx.SE(ix),FEx.Lower(ix),FEx.Upper(ix),FEx.pValue(ix), ...
    'VariableNames',{'Analysis','RemovedRID','Observations','Participants', ...
                     'TauTime_Beta','SE','CI_Lower','CI_Upper','PValue'});
InfluenceResults = [InfluenceResults; row];

fprintf('\nINFLUENCE SENSITIVITY\n');
disp(InfluenceResults);

%% Save outputs
writetable(ObsAudit,fullfile(outputDir,'Step13_ObservationResidual_Audit.csv'));
writetable(SubAudit,fullfile(outputDir,'Step13_ParticipantResidual_Audit.csv'));
writetable(InfluenceResults,fullfile(outputDir,'Step13_Influence_Sensitivity.csv'));

Summary = table(height(M),numel(unique(M.RID)),rMean,rMedian,rSD,rSkew,rKurt, ...
    nAbs3,rho,pRho,primaryBeta,primaryLo,primaryHi,primaryP, ...
    'VariableNames',{'Observations','Participants','ResidualMean','ResidualMedian', ...
    'ResidualSD','Skewness','Kurtosis','N_AbsStdResidual_GT3', ...
    'SpearmanAbsResidual_vs_Fitted','SpearmanP', ...
    'PrimaryTauTime_Beta','PrimaryTauTime_CI_Lower', ...
    'PrimaryTauTime_CI_Upper','PrimaryTauTime_P'});

writetable(Summary,fullfile(outputDir,'Step13_ResidualDiagnostics_Summary.csv'));

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');
