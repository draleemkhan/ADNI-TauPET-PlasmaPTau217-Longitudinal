%% ============================================================
% STEP 48B: FIXED-EFFECT MULTICOLLINEARITY / VIF DIAGNOSTICS
% =============================================================
%
% Purpose:
%   Assess multicollinearity in the fixed-effect structure of:
%
%   A. Primary model using RAW p-tau217
%   B. Sensitivity model using LOG p-tau217
%
%   VIF is calculated directly:
%
%       VIF_j = 1 / (1 - R_j^2)
%
%   where predictor j is regressed on all remaining predictors.
%
%   Particular interest:
%       Tau_z
%       PTau217_z / LogPTau217_z
%       TimeYears:Tau_z
%       TimeYears:PTau217_z / TimeYears:LogPTau217_z
%
% ============================================================

clear;
clc;

fprintf('\n============================================================\n');
fprintf(' STEP 48B: FIXED-EFFECT VIF DIAGNOSTICS\n');
fprintf('============================================================\n\n');

%% ------------------------------------------------------------
% Paths
% -------------------------------------------------------------

outputDir = ...
    'D:\output';

dataFile = fullfile(outputDir,...
    'ADNI_REVISED_COHORT_LOCKED_N143.mat');

assert(exist(dataFile,'file')==2,...
    'Locked analytic dataset not found.');

%% ------------------------------------------------------------
% Load analytic table
% -------------------------------------------------------------

tmp = load(dataFile);
fn = fieldnames(tmp);

T = [];

for k = 1:numel(fn)

    obj = tmp.(fn{k});

    if istable(obj)

        vars = string(obj.Properties.VariableNames);

        requiredCore = [ ...
            "RID","TimeYears","Age","Education",...
            "Tau_z","Centiloid_z",...
            "PTau217_Fujirebio"];

        if all(ismember(requiredCore,vars))
            T = obj;
            break;
        end
    end
end

assert(~isempty(T),...
    'Could not locate analytic table.');

fprintf('Dataset: %d observations\n',height(T));
fprintf('Participants: %d\n\n',numel(unique(T.RID)));

assert(height(T)==379,...
    'Expected 379 observations.');

assert(numel(unique(T.RID))==143,...
    'Expected 143 participants.');

fprintf('Dataset audit: PASS\n\n');

%% ------------------------------------------------------------
% Identify sex variable
% -------------------------------------------------------------

vars = string(T.Properties.VariableNames);

if ismember("SexModel",vars)

    sexRaw = string(T.SexModel);

elseif ismember("Sex_cat",vars)

    sexRaw = string(T.Sex_cat);

elseif ismember("Sex",vars)

    sexRaw = string(T.Sex);

else

    error('Could not identify sex variable.');

end

% Female = 1, Male = 0
SexFemale = double(strcmpi(sexRaw,'Female'));

fprintf('Sex coding: Male = 0, Female = 1\n');

%% ------------------------------------------------------------
% Construct log p-tau217 standardisation
% -------------------------------------------------------------

rawPTau = double(T.PTau217_Fujirebio);

assert(all(rawPTau > 0),...
    'Non-positive p-tau217 value encountered.');

logPTau = log(rawPTau);

% Standardise log p-tau217 using participant-level baseline
% distribution, matching Step 46B.

Tsort = sortrows(T,'RID');
[~,firstIdx] = unique(Tsort.RID,'stable');

baselineLogPTau = ...
    log(double(Tsort.PTau217_Fujirebio(firstIdx)));

muLog = mean(baselineLogPTau,'omitnan');
sdLog = std(baselineLogPTau,'omitnan');

LogPTau217_z = (logPTau - muLog) ./ sdLog;

fprintf('\nLOG p-tau217 standardisation\n');
fprintf('Participant-level baseline N: %d\n',...
    numel(baselineLogPTau));
fprintf('Mean log p-tau217: %.6f\n',muLog);
fprintf('SD log p-tau217:   %.6f\n\n',sdLog);

%% ------------------------------------------------------------
% Extract predictors
% -------------------------------------------------------------

TimeYears   = double(T.TimeYears);
Age         = double(T.Age);
Education   = double(T.Education);
Tau_z       = double(T.Tau_z);
Centiloid_z = double(T.Centiloid_z);

% Use existing standardized raw p-tau217 if available.
if ismember("PTau217_z",vars)

    PTau217_z = double(T.PTau217_z);

else

    % Fallback standardisation using participant-level baseline values
    baselineRaw = ...
        double(Tsort.PTau217_Fujirebio(firstIdx));

    muRaw = mean(baselineRaw,'omitnan');
    sdRaw = std(baselineRaw,'omitnan');

    PTau217_z = ...
        (double(T.PTau217_Fujirebio)-muRaw)./sdRaw;
end

%% ------------------------------------------------------------
% Construct interaction terms
% -------------------------------------------------------------

TauTime = TimeYears .* Tau_z;
CentiloidTime = TimeYears .* Centiloid_z;

RawPTauTime = TimeYears .* PTau217_z;
LogPTauTime = TimeYears .* LogPTau217_z;

%% ------------------------------------------------------------
% RAW p-tau217 fixed-effect design matrix
% -------------------------------------------------------------

Xraw = [ ...
    TimeYears,...
    Age,...
    Education,...
    Tau_z,...
    Centiloid_z,...
    PTau217_z,...
    SexFemale,...
    TauTime,...
    CentiloidTime,...
    RawPTauTime];

rawNames = { ...
    'TimeYears',...
    'Age',...
    'Education',...
    'Tau_z',...
    'Centiloid_z',...
    'PTau217_z',...
    'SexModel_Female',...
    'TimeYears_x_Tau_z',...
    'TimeYears_x_Centiloid_z',...
    'TimeYears_x_PTau217_z'};

%% ------------------------------------------------------------
% LOG p-tau217 fixed-effect design matrix
% -------------------------------------------------------------

Xlog = [ ...
    TimeYears,...
    Age,...
    Education,...
    Tau_z,...
    Centiloid_z,...
    LogPTau217_z,...
    SexFemale,...
    TauTime,...
    CentiloidTime,...
    LogPTauTime];

logNames = { ...
    'TimeYears',...
    'Age',...
    'Education',...
    'Tau_z',...
    'Centiloid_z',...
    'LogPTau217_z',...
    'SexModel_Female',...
    'TimeYears_x_Tau_z',...
    'TimeYears_x_Centiloid_z',...
    'TimeYears_x_LogPTau217_z'};

%% ------------------------------------------------------------
% Complete-case audit
% -------------------------------------------------------------

completeRaw = all(isfinite(Xraw),2);
completeLog = all(isfinite(Xlog),2);

fprintf('COMPLETE-CASE AUDIT\n');
fprintf('Raw model: %d / %d observations\n',...
    sum(completeRaw),height(T));
fprintf('Log model: %d / %d observations\n\n',...
    sum(completeLog),height(T));

assert(sum(completeRaw)==379,...
    'Raw model has incomplete predictor data.');

assert(sum(completeLog)==379,...
    'Log model has incomplete predictor data.');

%% ------------------------------------------------------------
% Calculate VIFs
% -------------------------------------------------------------

VIFraw = calculateVIF(Xraw(completeRaw,:));
VIFlog = calculateVIF(Xlog(completeLog,:));

RawVIF = table( ...
    string(rawNames(:)),...
    VIFraw(:),...
    'VariableNames',{'Term','VIF'});

LogVIF = table( ...
    string(logNames(:)),...
    VIFlog(:),...
    'VariableNames',{'Term','VIF'});

%% ------------------------------------------------------------
% Display full results
% -------------------------------------------------------------

fprintf('============================================================\n');
fprintf(' RAW p-tau217 MODEL VIFs\n');
fprintf('============================================================\n\n');

disp(RawVIF);

fprintf('Maximum VIF: %.4f\n\n',max(VIFraw));

fprintf('============================================================\n');
fprintf(' LOG p-tau217 MODEL VIFs\n');
fprintf('============================================================\n\n');

disp(LogVIF);

fprintf('Maximum VIF: %.4f\n\n',max(VIFlog));

%% ------------------------------------------------------------
% Key biomarker terms
% -------------------------------------------------------------

rawKeyIdx = [4 6 8 10];
logKeyIdx = [4 6 8 10];

KeyComparison = table( ...
    ["Tau PET main effect"; ...
     "p-tau217 main effect"; ...
     "Tau PET x Time"; ...
     "p-tau217 x Time"],...
    VIFraw(rawKeyIdx),...
    VIFlog(logKeyIdx),...
    VIFlog(logKeyIdx)-VIFraw(rawKeyIdx),...
    100.*(VIFlog(logKeyIdx)-VIFraw(rawKeyIdx))...
        ./VIFraw(rawKeyIdx),...
    'VariableNames',{ ...
    'Component',...
    'RawPTau_VIF',...
    'LogPTau_VIF',...
    'AbsoluteChange',...
    'PercentChange'});

fprintf('============================================================\n');
fprintf(' KEY RAW-vs-LOG VIF COMPARISON\n');
fprintf('============================================================\n\n');

disp(KeyComparison);

%% ------------------------------------------------------------
% Focus specifically on Tau PET x Time
% -------------------------------------------------------------

fprintf('\nTAU PET x TIME VIF\n');
fprintf('Raw p-tau217 model: %.4f\n',VIFraw(8));
fprintf('Log p-tau217 model: %.4f\n',VIFlog(8));
fprintf('Absolute change:     %.4f\n',...
    VIFlog(8)-VIFraw(8));

fprintf('Percent change:      %.2f%%\n',...
    100*(VIFlog(8)-VIFraw(8))/VIFraw(8));

%% ------------------------------------------------------------
% Condition indices as secondary global diagnostic
% -------------------------------------------------------------

% Standardize columns before singular-value decomposition.
Zraw = zscore(Xraw(completeRaw,:));
Zlog = zscore(Xlog(completeLog,:));

[~,Sraw,~] = svd(Zraw,0);
[~,Slog,~] = svd(Zlog,0);

sRaw = diag(Sraw);
sLog = diag(Slog);

condIndexRaw = max(sRaw)./sRaw;
condIndexLog = max(sLog)./sLog;

maxCondRaw = max(condIndexRaw);
maxCondLog = max(condIndexLog);

fprintf('\n============================================================\n');
fprintf(' GLOBAL CONDITION-INDEX DIAGNOSTIC\n');
fprintf('============================================================\n\n');

fprintf('Maximum condition index, RAW model: %.4f\n',...
    maxCondRaw);

fprintf('Maximum condition index, LOG model: %.4f\n',...
    maxCondLog);

%% ------------------------------------------------------------
% Save outputs
% -------------------------------------------------------------

writetable(RawVIF,...
    fullfile(outputDir,...
    'Step48B_RawPTau_FixedEffect_VIF.csv'));

writetable(LogVIF,...
    fullfile(outputDir,...
    'Step48B_LogPTau_FixedEffect_VIF.csv'));

writetable(KeyComparison,...
    fullfile(outputDir,...
    'Step48B_Key_VIF_Comparison.csv'));

ConditionSummary = table( ...
    ["Raw p-tau217";"Log p-tau217"],...
    [max(VIFraw);max(VIFlog)],...
    [maxCondRaw;maxCondLog],...
    'VariableNames',{ ...
    'Model',...
    'MaximumVIF',...
    'MaximumConditionIndex'});

writetable(ConditionSummary,...
    fullfile(outputDir,...
    'Step48B_ConditionIndex_Summary.csv'));

save(fullfile(outputDir,...
    'Step48B_VIF_Diagnostics.mat'),...
    'Xraw','Xlog',...
    'VIFraw','VIFlog',...
    'RawVIF','LogVIF',...
    'KeyComparison',...
    'condIndexRaw','condIndexLog',...
    'ConditionSummary');

%% ------------------------------------------------------------
% Complete
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' STEP 48B COMPLETE\n');
fprintf('============================================================\n');

fprintf('Observations: 379: PASS\n');
fprintf('Participants: 143: PASS\n');

fprintf('\nKEY DIAGNOSTIC\n');
fprintf('Tau PET x Time VIF, raw model: %.4f\n',VIFraw(8));
fprintf('Tau PET x Time VIF, log model: %.4f\n',VIFlog(8));

fprintf('\nMaximum VIF, raw model: %.4f\n',max(VIFraw));
fprintf('Maximum VIF, log model: %.4f\n',max(VIFlog));

fprintf('\nMaximum condition index, raw model: %.4f\n',...
    maxCondRaw);
fprintf('Maximum condition index, log model: %.4f\n',...
    maxCondLog);

fprintf('\nSaved outputs to:\n%s\n',outputDir);

fprintf('\nFiles created:\n');
fprintf('  Step48B_RawPTau_FixedEffect_VIF.csv\n');
fprintf('  Step48B_LogPTau_FixedEffect_VIF.csv\n');
fprintf('  Step48B_Key_VIF_Comparison.csv\n');
fprintf('  Step48B_ConditionIndex_Summary.csv\n');
fprintf('  Step48B_VIF_Diagnostics.mat\n');

fprintf('============================================================\n');

%% ============================================================
% LOCAL FUNCTION
% =============================================================

function vif = calculateVIF(X)

    p = size(X,2);
    vif = nan(p,1);

    for j = 1:p

        y = X(:,j);
        others = X(:,setdiff(1:p,j));

        % Include intercept in auxiliary regression
        Z = [ones(size(others,1),1), others];

        b = Z\y;
        yhat = Z*b;

        SSE = sum((y-yhat).^2);
        SST = sum((y-mean(y)).^2);

        R2 = 1 - SSE/SST;

        % Protect against numerical rounding above 1
        R2 = min(max(R2,0),1-eps);

        vif(j) = 1/(1-R2);

    end
end