%% Step12_RandomInterceptOnly_LME_Check.m
% Current locked ADNI cohort: N=143, 379 repeated CDR-SB observations.
% Purpose: fit random-intercept-only sensitivity model using the exact
% fixed-effects specification used in the locked primary model.

clear; clc; close all;

outputDir = 'D:\output';

% Prefer the exact model-ready dataset produced in Step08.
candidates = { ...
    fullfile(outputDir,'Step08_ModelReady_Longitudinal_Dataset.csv'), ...
    fullfile(outputDir,'Step08_ModelReady_Longitudinal.csv'), ...
    fullfile(outputDir,'Step08_Model_Ready_Longitudinal_Dataset.csv')};

modelFile = "";
for k=1:numel(candidates)
    if isfile(candidates{k}), modelFile=string(candidates{k}); break; end
end

% Fallback: identify a 379-row CSV containing required model variables.
if modelFile==""
    ff=dir(fullfile(outputDir,'*.csv'));
    needed=["RID","CDRSB","TimeYears","Tau_z","Centiloid_z","PTau217_z","Age","Education"];
    for k=1:numel(ff)
        try
            X=readtable(fullfile(ff(k).folder,ff(k).name),'VariableNamingRule','preserve');
            vn=string(X.Properties.VariableNames);
            if height(X)==379 && all(ismember(needed,vn))
                modelFile=string(fullfile(ff(k).folder,ff(k).name));
                break
            end
        catch
        end
    end
end
if modelFile=="", error('Could not locate the current 379-row Step08 model-ready dataset.'); end

T=readtable(modelFile,'VariableNamingRule','preserve');
fprintf('\n====================================================\n');
fprintf(' STEP 12 - RANDOM INTERCEPT ONLY LME CHECK\n');
fprintf('====================================================\n');
fprintf('Model-ready source: %s\n',modelFile);
fprintf('Observations: %d [expected 379]\n',height(T));
fprintf('Participants: %d [expected 143]\n',numel(unique(T.RID)));

% Reconstruct categorical variables only if needed.
if ~ismember('RIDcat',T.Properties.VariableNames)
    T.RIDcat=categorical(T.RID);
end
if ~ismember('SexModel',T.Properties.VariableNames)
    if ismember('Sex',T.Properties.VariableNames)
        T.SexModel=categorical(string(T.Sex));
    else
        error('Neither SexModel nor Sex found.');
    end
else
    T.SexModel=categorical(T.SexModel);
end

% Complete-case audit.
vars={'CDRSB','TimeYears','Tau_z','Centiloid_z','PTau217_z','Age','Education'};
ok=true(height(T),1);
for k=1:numel(vars)
    x=double(T.(vars{k}));
    ok=ok & isfinite(x);
end
ok=ok & ~isundefined(T.RIDcat) & ~isundefined(T.SexModel);
M=T(ok,:);

fprintf('Complete observations: %d\n',height(M));
fprintf('Complete participants: %d\n',numel(unique(M.RID)));

formula = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + Age + ' ...
           'SexModel + Education + Tau_z:TimeYears + Centiloid_z:TimeYears + ' ...
           'PTau217_z:TimeYears + (1|RIDcat)'];

fprintf('\nFORMULA\n%s\n',formula);

% REML: directly comparable to previously reported primary/sensitivity fits.
lme=fitlme(M,formula,'FitMethod','REML');

fprintf('\nMODEL FIT COMPLETE\n');
disp(lme);

FE=dataset2table(lme.Coefficients);
fprintf('\nFIXED EFFECTS\n');
disp(FE);

% Key Tau x Time
names=string(FE.Name);
idx=find(names=="TimeYears:Tau_z" | names=="Tau_z:TimeYears",1);
if isempty(idx), error('Tau x Time coefficient not found.'); end

fprintf('\n====================================================\n');
fprintf(' RANDOM-INTERCEPT-ONLY TAU x TIME\n');
fprintf('====================================================\n');
fprintf('Beta = %.6f\n',FE.Estimate(idx));
fprintf('SE = %.6f\n',FE.SE(idx));
fprintf('95%% CI = %.6f to %.6f\n',FE.Lower(idx),FE.Upper(idx));
fprintf('P = %.10g\n',FE.pValue(idx));

% Extract fit statistics robustly.
AIC=lme.ModelCriterion.AIC;
BIC=lme.ModelCriterion.BIC;
LL=lme.LogLikelihood;
fprintf('\nFIT STATISTICS\n');
fprintf('AIC = %.4f\n',AIC);
fprintf('BIC = %.4f\n',BIC);
fprintf('LogLikelihood = %.4f\n',LL);

% Save.
writetable(FE,fullfile(outputDir,'Step12_RandomInterceptOnly_FixedEffects.csv'));

S=table(height(M),numel(unique(M.RID)),AIC,BIC,LL, ...
    FE.Estimate(idx),FE.SE(idx),FE.Lower(idx),FE.Upper(idx),FE.pValue(idx), ...
    'VariableNames',{'Observations','Participants','AIC','BIC','LogLikelihood', ...
    'TauTime_Beta','TauTime_SE','TauTime_CI_Lower','TauTime_CI_Upper','TauTime_P'});
writetable(S,fullfile(outputDir,'Step12_RandomInterceptOnly_ModelSummary.csv'));

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');
