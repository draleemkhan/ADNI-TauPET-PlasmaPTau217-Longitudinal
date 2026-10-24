%% Step23_Participant_Cluster_Bootstrap_PrimaryTauTime.m
% ADNI revised cohort model
% Step 23: Participant-level cluster bootstrap of primary Tau x Time effect
%
% Resamples PARTICIPANTS (not rows) with replacement. All longitudinal
% observations for a selected participant travel together. Duplicate
% selections receive unique bootstrap cluster IDs.
%
% Primary model:
% CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + Age + SexModel +
%         Education + Tau_z:TimeYears + Centiloid_z:TimeYears +
%         PTau217_z:TimeYears
% Random effects: uncorrelated random intercept + random linear time slope.
%
% Bootstrap: 2000 replicates, REML.
%
% Outputs:
%   Step23_ClusterBootstrap_PrimaryTauTime.mat
%   Step23_ClusterBootstrap_Replicates.csv
%   Step23_ClusterBootstrap_Summary.csv

clear; clc;

%% Paths and settings
outDir = 'D:\output';
inFile = fullfile(outDir,'ADNI_REVISED_COHORT_LOCKED_N143.mat');
B = 2000;
rng(23092026,'twister'); % reproducible

assert(isfile(inFile),'Locked cohort file not found.');
S = load(inFile);
T = S.LockedLongitudinal;

fprintf('\nSTEP 23: PARTICIPANT-LEVEL CLUSTER BOOTSTRAP\n');
fprintf('============================================\n');
fprintf('Dataset: %d observations, %d participants\n',height(T),numel(unique(T.RID)));
fprintf('Requested bootstrap replicates: %d\n',B);

%% Locked primary formula
formula = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
           'Age + SexModel + Education + ' ...
           'Tau_z:TimeYears + Centiloid_z:TimeYears + PTau217_z:TimeYears + ' ...
           '(1|RIDcat) + (TimeYears-1|RIDcat)'];

%% Refit observed model for reference
obsModel = fitlme(T,formula,'FitMethod','REML','DummyVarCoding','reference');
obs = localGetCoef(obsModel,'Tau_z:TimeYears');

fprintf('\nObserved Tau x Time:\n');
fprintf('Beta = %.6f, SE = %.6f, 95%% CI %.6f to %.6f, P = %.6g\n', ...
    obs.Estimate,obs.SE,obs.Lower,obs.Upper,obs.pValue);

%% Prepare participant blocks
RIDlist = unique(T.RID,'stable');
N = numel(RIDlist);

Beta = nan(B,1);
SE = nan(B,1);
PValue = nan(B,1);
NObs = nan(B,1);
NOriginalRID = nan(B,1);
InterceptVar = nan(B,1);
SlopeVar = nan(B,1);
ResidualVar = nan(B,1);
NearBoundary = false(B,1);
FitSuccess = false(B,1);
FailureMessage = strings(B,1);

fprintf('\nRunning bootstrap...\n');

for b = 1:B
    draw = randi(N,N,1);

    % Build bootstrap table, assigning each sampled participant-copy
    % a unique grouping ID.
    parts = cell(N,1);
    for j = 1:N
        rid = RIDlist(draw(j));
        tmp = T(T.RID == rid,:);
        tmp.BootRID = repmat(categorical("B"+string(j)),height(tmp),1);
        parts{j} = tmp;
    end
    Tb = vertcat(parts{:});
    NObs(b) = height(Tb);
    NOriginalRID(b) = numel(unique(draw));

    % Replace grouping variable in formula with bootstrap cluster ID.
    fBoot = strrep(formula,'RIDcat','BootRID');

    try
        mdl = fitlme(Tb,fBoot,'FitMethod','REML','DummyVarCoding','reference');
        c = localGetCoef(mdl,'Tau_z:TimeYears');

        Beta(b) = c.Estimate;
        SE(b) = c.SE;
        PValue(b) = c.pValue;
        FitSuccess(b) = true;

        % Variance-component audit. covarianceParameters returns random
        % covariance parameter cells and residual SD/parameter depending
        % on MATLAB release/model structure. Use robust extraction from
        % model covariance outputs where possible.
        try
            [psi,mse] = covarianceParameters(mdl);
            ResidualVar(b) = mse;
            vals = [];
            for kk = 1:numel(psi)
                x = psi{kk};
                if isnumeric(x)
                    vals = [vals; x(:)]; %#ok<AGROW>
                end
            end
            if numel(vals) >= 2
                InterceptVar(b) = vals(1);
                SlopeVar(b) = vals(2);
                NearBoundary(b) = any(vals(1:2) < 1e-8) || mse < 1e-8;
            end
        catch
            % Do not fail an otherwise valid bootstrap fit solely because
            % variance-component extraction differs by MATLAB release.
        end

    catch ME
        FailureMessage(b) = string(ME.message);
    end

    if mod(b,100)==0 || b==1
        fprintf('  %4d/%d complete | successful fits: %d\n', ...
            b,B,sum(FitSuccess(1:b)));
    end
end

%% Summaries
ok = FitSuccess & isfinite(Beta);
nOK = sum(ok);
nFail = B-nOK;
successPct = 100*nOK/B;

if nOK < 0.90*B
    warning('Fewer than 90%% of bootstrap fits succeeded. Interpret cautiously.');
end

bootBeta = Beta(ok);
bootMean = mean(bootBeta);
bootMedian = median(bootBeta);
bootSD = std(bootBeta);
bootCI = prctile(bootBeta,[2.5 97.5]);
propPositive = mean(bootBeta > 0);
propNegative = mean(bootBeta < 0);
propNominalP = mean(PValue(ok) < 0.05);
bias = bootMean - obs.Estimate;

fprintf('\nBOOTSTRAP RESULTS\n');
fprintf('=================\n');
fprintf('Successful fits: %d/%d (%.1f%%)\n',nOK,B,successPct);
fprintf('Failed fits:     %d/%d (%.1f%%)\n',nFail,B,100*nFail/B);
fprintf('Observed beta:   %.6f\n',obs.Estimate);
fprintf('Bootstrap mean:  %.6f\n',bootMean);
fprintf('Bootstrap median:%.6f\n',bootMedian);
fprintf('Bootstrap SD:    %.6f\n',bootSD);
fprintf('Mean bias:       %.6f\n',bias);
fprintf('Percentile 95%% CI: %.6f to %.6f\n',bootCI(1),bootCI(2));
fprintf('Proportion beta > 0: %.4f (%.1f%%)\n',propPositive,100*propPositive);
fprintf('Proportion beta < 0: %.4f (%.1f%%)\n',propNegative,100*propNegative);
fprintf('Proportion nominal model P < .05: %.4f (%.1f%%)\n', ...
    propNominalP,100*propNominalP);
fprintf('Near-boundary variance fits detected: %d\n',sum(NearBoundary(ok)));

%% Replicate table
Replicates = table((1:B)',FitSuccess,Beta,SE,PValue,NObs,NOriginalRID, ...
    InterceptVar,SlopeVar,ResidualVar,NearBoundary,FailureMessage, ...
    'VariableNames',{'Replicate','FitSuccess','TauTime_Beta','TauTime_SE', ...
    'TauTime_P','Observations','UniqueOriginalParticipants','InterceptComponent', ...
    'SlopeComponent','ResidualComponent','NearBoundary','FailureMessage'});

Summary = table( ...
    B,nOK,nFail,successPct,obs.Estimate,obs.SE,obs.Lower,obs.Upper,obs.pValue, ...
    bootMean,bootMedian,bootSD,bias,bootCI(1),bootCI(2), ...
    propPositive,propNegative,propNominalP,sum(NearBoundary(ok)), ...
    'VariableNames',{'RequestedReplicates','SuccessfulFits','FailedFits', ...
    'SuccessPercent','ObservedBeta','ObservedSE','ObservedCI_Lower', ...
    'ObservedCI_Upper','ObservedP','BootstrapMean','BootstrapMedian', ...
    'BootstrapSD','MeanBias','BootstrapCI_Lower','BootstrapCI_Upper', ...
    'ProportionPositive','ProportionNegative','ProportionNominalP_LT_05', ...
    'NearBoundaryFits'});

disp(Summary);

%% Save
matFile = fullfile(outDir,'Step23_ClusterBootstrap_PrimaryTauTime.mat');
save(matFile,'Replicates','Summary','Beta','SE','PValue','FitSuccess', ...
    'NearBoundary','obs','formula','B','-v7.3');

writetable(Replicates,fullfile(outDir,'Step23_ClusterBootstrap_Replicates.csv'));
writetable(Summary,fullfile(outDir,'Step23_ClusterBootstrap_Summary.csv'));

fprintf('\nSTEP 23 COMPLETE\n');
fprintf('Saved MAT and CSV outputs to:\n%s\n',outDir);
fprintf('\nPlease paste the full final MATLAB output for interpretation.\n');

%% Local helper
function out = localGetCoef(mdl,coefName)
    C = mdl.Coefficients;
    idx = strcmp(C.Name,coefName);
    if ~any(idx) && contains(coefName,':')
        z = split(string(coefName),':');
        alt = char(z(2)+":"+z(1));
        idx = strcmp(C.Name,alt);
    end
    if ~any(idx)
        error('Coefficient "%s" not found.',coefName);
    end
    r = C(find(idx,1),:);
    out = struct('Estimate',r.Estimate,'SE',r.SE,'Lower',r.Lower, ...
                 'Upper',r.Upper,'pValue',r.pValue);
end
