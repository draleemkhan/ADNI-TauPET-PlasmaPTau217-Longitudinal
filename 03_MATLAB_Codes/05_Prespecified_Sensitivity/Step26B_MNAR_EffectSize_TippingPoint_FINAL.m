%% Step26B_MNAR_EffectSize_TippingPoint_FINAL.m
% ADNI revised cohort
% FINAL Step 26B: Controlled MNAR effect-size tipping-point sensitivity
%
% PURPOSE
% Quantify how different the Tau_z x TimeYears association among the
% 67 non-retained participants would have to be to eliminate the overall
% point estimate across all 210 baseline-eligible participants.
%
% This analysis deliberately DOES NOT:
%   - fabricate missing longitudinal CDR-SB trajectories;
%   - assign P values to hypothetical missing-participant effects;
%   - claim to identify the true effect among non-retained participants.
%
% FRAMEWORK
% beta_O = IPCW-adjusted Tau x Time estimate among retained participants.
% beta_M = assumed Tau x Time effect among non-retained participants.
%
% beta_overall = pi_O*beta_O + pi_M*beta_M
%
% The sensitivity parameter can also be written as:
%   k = beta_M / beta_O
%
% where k=1 means the same effect, k=0 means no effect, k=-1 means an
% equal-magnitude association in the opposite direction.
%
% The tipping point is the beta_M (and k) for which beta_overall = 0.

clear; clc;

%% ---------------- PATHS ----------------
outDir = 'D:\output';
step26AFile = fullfile(outDir,'Step26A_IPCW_Retention_Sensitivity.mat');

assert(isfile(step26AFile),'Cannot find Step26A MAT file: %s',step26AFile);
S = load(step26AFile);

%% ---------------- COUNTS ----------------
Ntotal = 210;
Nobs   = 143;
Nmiss  = 67;

piO = Nobs/Ntotal;
piM = Nmiss/Ntotal;

%% ---------------- IPCW ANCHOR ----------------
R = S.Result;
names = string(R.Model);
idx = find(contains(lower(names),'ipcw untrimmed'),1);
assert(~isempty(idx),'Could not find IPCW untrimmed estimate in Step26A Result.');

betaO = R.TauTimeBeta(idx);

fprintf('\n============================================================\n');
fprintf(' STEP 26B FINAL: MNAR EFFECT-SIZE TIPPING POINT\n');
fprintf('============================================================\n');
fprintf('Baseline eligible N: %d\n',Ntotal);
fprintf('Retained N: %d (%.1f%%)\n',Nobs,100*piO);
fprintf('Non-retained N: %d (%.1f%%)\n',Nmiss,100*piM);
fprintf('IPCW-adjusted retained-group Tau x Time beta: %.6f\n',betaO);

%% ---------------- EXACT POINT-ESTIMATE TIPPING POINT ----------------
% Solve:
% 0 = piO*betaO + piM*betaM_tip

betaM_tip = -(piO/piM)*betaO;
k_tip = betaM_tip/betaO;

fprintf('\nExact point-estimate tipping point:\n');
fprintf('Assumed missing-group beta required for overall beta = 0: %.6f\n',betaM_tip);
fprintf('Relative to observed IPCW beta (k): %.6f\n',k_tip);
fprintf('Interpretation: %.2f times the observed magnitude in the opposite direction.\n',abs(k_tip));

%% ---------------- INTERPRETABLE SCENARIOS ----------------
Scenario = [
    "Same effect as retained (k = 1.0)"
    "75% of retained effect (k = 0.75)"
    "50% of retained effect (k = 0.50)"
    "25% of retained effect (k = 0.25)"
    "No tau effect (k = 0)"
    "Half magnitude, opposite direction (k = -0.50)"
    "Equal magnitude, opposite direction (k = -1.0)"
    "Twice magnitude, opposite direction (k = -2.0)"
    "Exact point-estimate tipping point"
    ];

k = [1; .75; .50; .25; 0; -.50; -1; -2; k_tip];
betaM = k.*betaO;
betaOverall = piO.*betaO + piM.*betaM;
percentOfAnchor = 100.*betaOverall./betaO;

Scenarios = table(Scenario,k,betaM,betaOverall,percentOfAnchor, ...
    'VariableNames',{'Scenario','RelativeEffect_k','AssumedMissingBeta', ...
                     'OverallMixtureBeta','OverallPercentOfIPCWAnchor'});

fprintf('\n============================================================\n');
fprintf(' INTERPRETABLE MNAR EFFECT-SIZE SCENARIOS\n');
fprintf('============================================================\n');
disp(Scenarios);

%% ---------------- FINE EFFECT-DEPARTURE GRID ----------------
% k from +1.5 to -3.0. This is saved for complete reproducibility and
% possible supplementary plotting.

kGrid = (1.50:-0.01:-3.00)';
betaMGrid = kGrid.*betaO;
betaOverallGrid = piO.*betaO + piM.*betaMGrid;
percentAnchorGrid = 100.*betaOverallGrid./betaO;

Grid = table(kGrid,betaMGrid,betaOverallGrid,percentAnchorGrid, ...
    'VariableNames',{'RelativeEffect_k','AssumedMissingBeta', ...
                     'OverallMixtureBeta','OverallPercentOfIPCWAnchor'});

[~,nearTip] = min(abs(Grid.OverallMixtureBeta));

fprintf('\n============================================================\n');
fprintf(' GRID CHECK NEAREST TO ZERO\n');
fprintf('============================================================\n');
disp(Grid(nearTip,:));

%% ---------------- SUMMARY ----------------
Summary = table( ...
    Ntotal,Nobs,Nmiss,piO,piM,betaO,betaM_tip,k_tip, ...
    'VariableNames',{'BaselineEligibleN','RetainedN','NonRetainedN', ...
    'RetainedProportion','NonRetainedProportion','IPCW_AnchorBeta', ...
    'MissingBeta_PointEstimateTip','RelativeEffect_k_Tip'});

fprintf('\n============================================================\n');
fprintf(' FINAL SUMMARY\n');
fprintf('============================================================\n');
disp(Summary);

%% ---------------- SAVE ----------------
writetable(Scenarios,fullfile(outDir,'Step26B_FINAL_MNAR_EffectSize_Scenarios.csv'));
writetable(Grid,fullfile(outDir,'Step26B_FINAL_MNAR_EffectSize_Grid.csv'));
writetable(Summary,fullfile(outDir,'Step26B_FINAL_MNAR_EffectSize_Summary.csv'));

save(fullfile(outDir,'Step26B_FINAL_MNAR_EffectSize_TippingPoint.mat'), ...
    'Scenarios','Grid','Summary','Ntotal','Nobs','Nmiss','piO','piM', ...
    'betaO','betaM_tip','k_tip');

fprintf('\n============================================================\n');
fprintf(' STEP 26B FINAL COMPLETE\n');
fprintf('============================================================\n');
fprintf('Saved outputs to:\n%s\n',outDir);
fprintf('\nPlease paste:\n');
fprintf('1) Exact point-estimate tipping point\n');
fprintf('2) INTERPRETABLE MNAR EFFECT-SIZE SCENARIOS\n');
fprintf('3) FINAL SUMMARY\n');
