%% Step30B_Bootstrap_Incremental_ModelFit_R2Corrected.m
% ADNI revised cohort model
% Step 30B: Corrected participant-level bootstrap stability of incremental model fit
%
% Purpose:
%   Assess sampling stability of the incremental IN-SAMPLE contribution
%   of tau PET beyond demographics, Centiloid, and plasma p-tau217.
%
% IMPORTANT:
%   This is NOT external validation and is NOT labelled validated prediction.
%
% M2: demographics + Centiloid + p-tau217 (main effects + time interactions)
% M3: M2 + tau PET main effect + Tau x Time
%
% For each participant-level bootstrap sample:
%   1) fit M2 and M3 using REML -> marginal/conditional R2
%   2) calculate Delta marginal R2 = M3 - M2
%   3) fit M2 and M3 using ML -> logLik, AIC, BIC
%   4) calculate LR statistic, Delta AIC and Delta BIC
%
% Bootstrap resamples participants with replacement, preserving all repeated
% observations. Duplicate participants receive unique bootstrap cluster IDs.
%
% R2 estimator harmonized exactly to Step 19.
% Outputs saved to D:\output.

clear; clc;

%% Paths/settings
outDir = 'D:\output';
inFile = fullfile(outDir,'ADNI_REVISED_COHORT_LOCKED_N143.mat');
B = 2000;
rng(25092026,'twister');

assert(isfile(inFile),'Locked cohort file not found.');
S = load(inFile);
T = S.LockedLongitudinal;

fprintf('\nSTEP 25: BOOTSTRAP STABILITY OF INCREMENTAL MODEL FIT\n');
fprintf('====================================================\n');
fprintf('Dataset: %d observations, %d participants\n',height(T),numel(unique(T.RID)));
fprintf('Bootstrap replicates: %d\n',B);

%% Exact M2/M3 model definitions
fixedM2 = ['CDRSB ~ TimeYears + Centiloid_z + PTau217_z + Age + SexModel + Education + ' ...
           'Centiloid_z:TimeYears + PTau217_z:TimeYears'];

fixedM3 = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + Age + SexModel + Education + ' ...
           'Tau_z:TimeYears + Centiloid_z:TimeYears + PTau217_z:TimeYears'];

randPart = ' + (1|RIDcat) + (TimeYears-1|RIDcat)';
formM2 = [fixedM2 randPart];
formM3 = [fixedM3 randPart];

%% Original-sample fits
fprintf('\nFitting original-sample models...\n');

M2_REML = fitlme(T,formM2,'FitMethod','REML','DummyVarCoding','reference');
M3_REML = fitlme(T,formM3,'FitMethod','REML','DummyVarCoding','reference');

[r2m2,r2c2] = nakagawaR2(M2_REML);
[r2m3,r2c3] = nakagawaR2(M3_REML);

obsDeltaR2 = r2m3-r2m2;

M2_ML = fitlme(T,formM2,'FitMethod','ML','DummyVarCoding','reference');
M3_ML = fitlme(T,formM3,'FitMethod','ML','DummyVarCoding','reference');

obsLL2 = M2_ML.LogLikelihood;
obsLL3 = M3_ML.LogLikelihood;
obsLR = 2*(obsLL3-obsLL2);
obsDeltaAIC = M3_ML.ModelCriterion.AIC-M2_ML.ModelCriterion.AIC;
obsDeltaBIC = M3_ML.ModelCriterion.BIC-M2_ML.ModelCriterion.BIC;

fprintf('\nORIGINAL SAMPLE\n');
fprintf('M2 marginal R2 = %.6f | conditional R2 = %.6f\n',r2m2,r2c2);
fprintf('M3 marginal R2 = %.6f | conditional R2 = %.6f\n',r2m3,r2c3);
fprintf('Delta marginal R2 (M3-M2) = %.6f\n',obsDeltaR2);
fprintf('ML LR statistic = %.6f (df=2)\n',obsLR);
fprintf('Delta AIC (M3-M2) = %.6f\n',obsDeltaAIC);
fprintf('Delta BIC (M3-M2) = %.6f\n',obsDeltaBIC);

%% Bootstrap
RIDlist = unique(T.RID,'stable');
N = numel(RIDlist);

R2M2 = nan(B,1); R2M3 = nan(B,1); DeltaR2 = nan(B,1);
CondR2M2 = nan(B,1); CondR2M3 = nan(B,1);
LL2 = nan(B,1); LL3 = nan(B,1); LRStat = nan(B,1);
DeltaAIC = nan(B,1); DeltaBIC = nan(B,1);
FitSuccess = false(B,1);
FailureMessage = strings(B,1);
NObs = nan(B,1); UniqueOriginalParticipants = nan(B,1);

fprintf('\nRunning participant-level bootstrap...\n');

for b=1:B
    draw = randi(N,N,1);
    parts = cell(N,1);

    for j=1:N
        rid = RIDlist(draw(j));
        tmp = T(T.RID==rid,:);
        tmp.BootRID = repmat(categorical("B"+string(j)),height(tmp),1);
        parts{j}=tmp;
    end
    Tb = vertcat(parts{:});
    NObs(b)=height(Tb);
    UniqueOriginalParticipants(b)=numel(unique(draw));

    f2 = strrep(formM2,'RIDcat','BootRID');
    f3 = strrep(formM3,'RIDcat','BootRID');

    try
        % REML for descriptive marginal/conditional R2
        m2r = fitlme(Tb,f2,'FitMethod','REML','DummyVarCoding','reference');
        m3r = fitlme(Tb,f3,'FitMethod','REML','DummyVarCoding','reference');

        [R2M2(b),CondR2M2(b)] = nakagawaR2(m2r);
        [R2M3(b),CondR2M3(b)] = nakagawaR2(m3r);
        DeltaR2(b)=R2M3(b)-R2M2(b);

        % ML for valid nested fixed-effect model comparison
        m2m = fitlme(Tb,f2,'FitMethod','ML','DummyVarCoding','reference');
        m3m = fitlme(Tb,f3,'FitMethod','ML','DummyVarCoding','reference');

        LL2(b)=m2m.LogLikelihood;
        LL3(b)=m3m.LogLikelihood;
        LRStat(b)=2*(LL3(b)-LL2(b));
        DeltaAIC(b)=m3m.ModelCriterion.AIC-m2m.ModelCriterion.AIC;
        DeltaBIC(b)=m3m.ModelCriterion.BIC-m2m.ModelCriterion.BIC;

        FitSuccess(b)=all(isfinite([R2M2(b),R2M3(b),DeltaR2(b), ...
            LL2(b),LL3(b),LRStat(b),DeltaAIC(b),DeltaBIC(b)]));
    catch ME
        FailureMessage(b)=string(ME.message);
    end

    if mod(b,100)==0 || b==1
        fprintf('  %4d/%d complete | successful: %d\n', ...
            b,B,sum(FitSuccess(1:b)));
    end
end

%% Summaries
ok=FitSuccess;
nOK=sum(ok); nFail=B-nOK;

dR=DeltaR2(ok);
lr=LRStat(ok);
da=DeltaAIC(ok);
db=DeltaBIC(ok);

ciDR=prctile(dR,[2.5 97.5]);
ciLR=prctile(lr,[2.5 97.5]);
ciDA=prctile(da,[2.5 97.5]);
ciDB=prctile(db,[2.5 97.5]);

propDeltaPositive=mean(dR>0);
propAICImproved=mean(da<0);
propBICImproved=mean(db<0);
propLR_p_lt05=mean(1-chi2cdf(lr,2)<0.05);

fprintf('\nBOOTSTRAP RESULTS\n');
fprintf('=================\n');
fprintf('Successful fits: %d/%d (%.1f%%)\n',nOK,B,100*nOK/B);
fprintf('Failed fits: %d/%d (%.1f%%)\n',nFail,B,100*nFail/B);

fprintf('\nDelta marginal R2 (M3-M2):\n');
fprintf('Observed: %.6f\n',obsDeltaR2);
fprintf('Bootstrap mean: %.6f\n',mean(dR));
fprintf('Bootstrap median: %.6f\n',median(dR));
fprintf('Bootstrap SD: %.6f\n',std(dR));
fprintf('Percentile 95%% interval: %.6f to %.6f\n',ciDR(1),ciDR(2));
fprintf('Proportion Delta R2 > 0: %.4f (%.1f%%)\n',propDeltaPositive,100*propDeltaPositive);

fprintf('\nML nested model-fit improvement:\n');
fprintf('Observed LR statistic: %.6f\n',obsLR);
fprintf('Bootstrap median LR: %.6f\n',median(lr));
fprintf('Bootstrap LR 95%% interval: %.6f to %.6f\n',ciLR(1),ciLR(2));
fprintf('Proportion nominal LR P<.05 (df=2): %.4f (%.1f%%)\n', ...
    propLR_p_lt05,100*propLR_p_lt05);

fprintf('\nDelta AIC (M3-M2; negative favours M3):\n');
fprintf('Observed: %.6f | bootstrap median: %.6f | 95%% interval %.6f to %.6f\n', ...
    obsDeltaAIC,median(da),ciDA(1),ciDA(2));
fprintf('Proportion Delta AIC < 0: %.4f (%.1f%%)\n',propAICImproved,100*propAICImproved);

fprintf('\nDelta BIC (M3-M2; negative favours M3):\n');
fprintf('Observed: %.6f | bootstrap median: %.6f | 95%% interval %.6f to %.6f\n', ...
    obsDeltaBIC,median(db),ciDB(1),ciDB(2));
fprintf('Proportion Delta BIC < 0: %.4f (%.1f%%)\n',propBICImproved,100*propBICImproved);

%% Save
Replicates=table((1:B)',FitSuccess,R2M2,R2M3,DeltaR2,CondR2M2,CondR2M3, ...
    LL2,LL3,LRStat,DeltaAIC,DeltaBIC,NObs,UniqueOriginalParticipants,FailureMessage, ...
    'VariableNames',{'Replicate','FitSuccess','M2_MarginalR2','M3_MarginalR2', ...
    'DeltaMarginalR2','M2_ConditionalR2','M3_ConditionalR2','M2_LogLik','M3_LogLik', ...
    'LRStat','DeltaAIC','DeltaBIC','Observations','UniqueOriginalParticipants', ...
    'FailureMessage'});

Summary=table(B,nOK,nFail,r2m2,r2m3,obsDeltaR2,mean(dR),median(dR),std(dR), ...
    ciDR(1),ciDR(2),propDeltaPositive,obsLR,median(lr),ciLR(1),ciLR(2), ...
    propLR_p_lt05,obsDeltaAIC,median(da),ciDA(1),ciDA(2),propAICImproved, ...
    obsDeltaBIC,median(db),ciDB(1),ciDB(2),propBICImproved, ...
    'VariableNames',{'RequestedReplicates','SuccessfulFits','FailedFits', ...
    'ObservedM2MarginalR2','ObservedM3MarginalR2','ObservedDeltaMarginalR2', ...
    'BootstrapMeanDeltaR2','BootstrapMedianDeltaR2','BootstrapSDDeltaR2', ...
    'BootstrapDeltaR2_CI_Lower','BootstrapDeltaR2_CI_Upper','ProportionDeltaR2Positive', ...
    'ObservedLR','BootstrapMedianLR','BootstrapLR_CI_Lower','BootstrapLR_CI_Upper', ...
    'ProportionLR_P_LT_05','ObservedDeltaAIC','BootstrapMedianDeltaAIC', ...
    'BootstrapDeltaAIC_CI_Lower','BootstrapDeltaAIC_CI_Upper','ProportionDeltaAIC_LT_0', ...
    'ObservedDeltaBIC','BootstrapMedianDeltaBIC','BootstrapDeltaBIC_CI_Lower', ...
    'BootstrapDeltaBIC_CI_Upper','ProportionDeltaBIC_LT_0'});

disp(Summary);

save(fullfile(outDir,'Step30B_Bootstrap_Incremental_ModelFit_R2Corrected.mat'), ...
    'Replicates','Summary','M2_REML','M3_REML','M2_ML','M3_ML', ...
    'formM2','formM3','B','-v7.3');
writetable(Replicates,fullfile(outDir,'Step30B_Bootstrap_Incremental_ModelFit_R2Corrected_Replicates.csv'));
writetable(Summary,fullfile(outDir,'Step30B_Bootstrap_Incremental_ModelFit_R2Corrected_Summary.csv'));

fprintf('\nSTEP 30B COMPLETE\n');
fprintf('Saved MAT and CSV outputs to:\n%s\n',outDir);
fprintf('\nPlease paste the ORIGINAL SAMPLE and BOOTSTRAP RESULTS sections.\n');

%% Local function: Nakagawa-style marginal/conditional R2
function [R2m,R2c]=nakagawaR2(mdl)
    % Same variance-component R2 definition as Step 19, adapted to the
    % equivalent diagonal random-effects syntax used in this bootstrap:
    %   (1|cluster) + (TimeYears-1|cluster)
    %
    % In MATLAB this syntax can return the random intercept and random slope
    % as TWO separate 1x1 covariance components rather than one 2x2 matrix.

    % Fixed-effect variance
    yFixed = fitted(mdl,'Conditional',false);
    varFixed = var(yFixed,1);

    % Extract model-estimated random-effect covariance component(s)
    [psi,~,~] = covarianceParameters(mdl);

    covIR = 0;  % diagonal/uncorrelated structure by construction

    if iscell(psi)
        if numel(psi) >= 2 && isscalar(psi{1}) && isscalar(psi{2})
            % Separate random-intercept and random-slope terms
            varRI = double(psi{1});
            varRS = double(psi{2});
        elseif numel(psi) >= 1 && size(psi{1},1) >= 2 && size(psi{1},2) >= 2
            % Single 2x2 covariance matrix
            G = double(psi{1});
            varRI = G(1,1);
            varRS = G(2,2);
            covIR = G(1,2);
        else
            error('Could not identify random-intercept and random-slope covariance components.');
        end
    else
        G = double(psi);
        if size(G,1) >= 2 && size(G,2) >= 2
            varRI = G(1,1);
            varRS = G(2,2);
            covIR = G(1,2);
        elseif numel(G) >= 2
            varRI = G(1);
            varRS = G(2);
        else
            error('Could not identify random-intercept and random-slope covariance components.');
        end
    end

    % Observation-specific random-effect variance z_i G z_i'
    t = double(mdl.Variables.TimeYears);
    randVarObs = varRI + 2*t*covIR + (t.^2)*varRS;
    varRandom = mean(randVarObs);

    % Residual variance
    varResidual = mdl.MSE;

    denom = varFixed + varRandom + varResidual;
    R2m = varFixed / denom;
    R2c = (varFixed + varRandom) / denom;
end
