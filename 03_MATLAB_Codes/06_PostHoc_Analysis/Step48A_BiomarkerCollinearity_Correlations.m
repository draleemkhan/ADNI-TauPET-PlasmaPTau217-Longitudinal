%% ============================================================
% STEP 48A: PARTICIPANT-LEVEL BIOMARKER COLLINEARITY AUDIT
% =============================================================
%
% Purpose:
%   Examine participant-level baseline correlations among:
%       1. Tau PET
%       2. Centiloid
%       3. Raw plasma p-tau217
%       4. Log plasma p-tau217
%
%   Particular interest:
%       Tau PET vs raw p-tau217
%       Tau PET vs log p-tau217
%
% Outputs:
%   Pearson correlation matrix
%   Spearman correlation matrix
%   Pairwise correlation table
%   CSV files saved to D:\output
%
% =============================================================

clear;
clc;

fprintf('\n============================================================\n');
fprintf(' STEP 48A: PARTICIPANT-LEVEL BIOMARKER CORRELATIONS\n');
fprintf('============================================================\n\n');

%% ------------------------------------------------------------
% Paths
% -------------------------------------------------------------

outputDir = ...
    'D:\output';

if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

%% ------------------------------------------------------------
% Locate analytic dataset
% -------------------------------------------------------------

candidateFiles = { ...
    fullfile(outputDir,'Step39A_FinalPrimaryLME_ForFigure2.mat'), ...
    fullfile(outputDir,'Step47_LevelVsSlope_Models.mat'), ...
    fullfile(outputDir,'Step46B_LogPTau217_Sensitivity_Models.mat')};

S = [];
loadedFile = '';

for k = 1:numel(candidateFiles)
    if exist(candidateFiles{k},'file')
        tmp = load(candidateFiles{k});

        f = fieldnames(tmp);

        for j = 1:numel(f)
            obj = tmp.(f{j});

            if istable(obj)
                vars = string(obj.Properties.VariableNames);

                if all(ismember( ...
                        ["RID","TauMetaSUVR", ...
                         "PTau217_Fujirebio"],vars))

                    S = obj;
                    loadedFile = candidateFiles{k};
                    break;
                end
            end
        end

        if ~isempty(S)
            break;
        end
    end
end

%% ------------------------------------------------------------
% Fallback: inspect MAT files in D:\output
% -------------------------------------------------------------

if isempty(S)

    matFiles = dir(fullfile(outputDir,'*.mat'));

    for k = 1:numel(matFiles)

        fileName = fullfile(matFiles(k).folder,...
                            matFiles(k).name);

        tmp = load(fileName);
        f = fieldnames(tmp);

        for j = 1:numel(f)

            obj = tmp.(f{j});

            if istable(obj)

                vars = string(obj.Properties.VariableNames);

                if all(ismember( ...
                        ["RID","TauMetaSUVR", ...
                         "PTau217_Fujirebio"],vars))

                    S = obj;
                    loadedFile = fileName;
                    break;
                end
            end
        end

        if ~isempty(S)
            break;
        end
    end
end

assert(~isempty(S), ...
    ['Could not locate the analytic table. ' ...
     'Load the final analytic dataset manually.']);

fprintf('Analytic table source:\n%s\n\n',loadedFile);

fprintf('Observation-level dataset: %d rows\n',height(S));
fprintf('Unique participants: %d\n\n',numel(unique(S.RID)));

%% ------------------------------------------------------------
% Identify Centiloid variable
% -------------------------------------------------------------

vars = string(S.Properties.VariableNames);

centiloidCandidates = [ ...
    "Centiloid", ...
    "CentiloidRaw", ...
    "Centiloid_Value", ...
    "CENTILOIDS", ...
    "Centiloid_z"];

centVar = "";

for k = 1:numel(centiloidCandidates)
    if ismember(centiloidCandidates(k),vars)
        centVar = centiloidCandidates(k);
        break;
    end
end

assert(centVar ~= "", ...
    'Could not identify a Centiloid variable.');

fprintf('Centiloid variable used: %s\n',centVar);

%% ------------------------------------------------------------
% Create one row per participant
% -------------------------------------------------------------

% Sort by participant ID. The biomarkers are baseline/fixed
% participant-level predictors in the analytic dataset.
S = sortrows(S,'RID');

[~,firstIdx] = unique(S.RID,'stable');

B = S(firstIdx,:);

fprintf('Participant-level rows retained: %d\n',height(B));

assert(height(B)==143, ...
    'Expected 143 participants.');

fprintf('Participant-count audit: PASS\n\n');

%% ------------------------------------------------------------
% Extract biomarkers
% -------------------------------------------------------------

TauPET = double(B.TauMetaSUVR);
PTauRaw = double(B.PTau217_Fujirebio);
Centiloid = double(B.(char(centVar)));

% If only standardized Centiloid is available, correlation is
% unchanged by linear standardisation.
PTauLog = log(PTauRaw);

assert(all(PTauRaw > 0), ...
    'p-tau217 contains non-positive values; log invalid.');

%% ------------------------------------------------------------
% Missing-value audit
% -------------------------------------------------------------

fprintf('MISSING-VALUE AUDIT\n');

fprintf('Tau PET:        %d missing\n',sum(isnan(TauPET)));
fprintf('Centiloid:      %d missing\n',sum(isnan(Centiloid)));
fprintf('Raw p-tau217:   %d missing\n',sum(isnan(PTauRaw)));
fprintf('Log p-tau217:   %d missing\n\n',sum(isnan(PTauLog)));

complete = ...
    ~isnan(TauPET) & ...
    ~isnan(Centiloid) & ...
    ~isnan(PTauRaw) & ...
    ~isnan(PTauLog);

fprintf('Complete participants for correlation analysis: %d\n\n',...
    sum(complete));

X = [ ...
    TauPET(complete), ...
    Centiloid(complete), ...
    PTauRaw(complete), ...
    PTauLog(complete)];

biomarkerNames = { ...
    'Tau PET', ...
    'Centiloid', ...
    'Raw p-tau217', ...
    'Log p-tau217'};

%% ------------------------------------------------------------
% Pearson correlations
% -------------------------------------------------------------

[Rpearson,Ppearson] = corr(X,...
    'Type','Pearson',...
    'Rows','complete');

fprintf('============================================================\n');
fprintf(' PEARSON CORRELATION MATRIX\n');
fprintf('============================================================\n\n');

PearsonTable = array2table(Rpearson,...
    'VariableNames', ...
    {'TauPET','Centiloid','RawPTau217','LogPTau217'},...
    'RowNames', ...
    {'TauPET','Centiloid','RawPTau217','LogPTau217'});

disp(PearsonTable);

%% ------------------------------------------------------------
% Spearman correlations
% -------------------------------------------------------------

[Rspearman,Pspearman] = corr(X,...
    'Type','Spearman',...
    'Rows','complete');

fprintf('\n============================================================\n');
fprintf(' SPEARMAN CORRELATION MATRIX\n');
fprintf('============================================================\n\n');

SpearmanTable = array2table(Rspearman,...
    'VariableNames',...
    {'TauPET','Centiloid','RawPTau217','LogPTau217'},...
    'RowNames',...
    {'TauPET','Centiloid','RawPTau217','LogPTau217'});

disp(SpearmanTable);

%% ------------------------------------------------------------
% Pairwise results table
% -------------------------------------------------------------

pairA = {};
pairB = {};
pearsonR = [];
pearsonP = [];
spearmanRho = [];
spearmanP = [];

row = 0;

for i = 1:size(X,2)-1
    for j = i+1:size(X,2)

        row = row + 1;

        pairA{row,1} = biomarkerNames{i};
        pairB{row,1} = biomarkerNames{j};

        pearsonR(row,1) = Rpearson(i,j);
        pearsonP(row,1) = Ppearson(i,j);

        spearmanRho(row,1) = Rspearman(i,j);
        spearmanP(row,1) = Pspearman(i,j);
    end
end

PairwiseResults = table( ...
    pairA,pairB,...
    pearsonR,pearsonP,...
    spearmanRho,spearmanP,...
    'VariableNames',{ ...
    'Biomarker1',...
    'Biomarker2',...
    'PearsonR',...
    'PearsonP',...
    'SpearmanRho',...
    'SpearmanP'});

fprintf('\n============================================================\n');
fprintf(' PAIRWISE BIOMARKER CORRELATIONS\n');
fprintf('============================================================\n\n');

disp(PairwiseResults);

%% ------------------------------------------------------------
% Key raw-vs-log comparison
% -------------------------------------------------------------

rTauRaw = Rpearson(1,3);
rTauLog = Rpearson(1,4);

rhoTauRaw = Rspearman(1,3);
rhoTauLog = Rspearman(1,4);

fprintf('\n============================================================\n');
fprintf(' KEY TAU PET / p-tau217 COMPARISON\n');
fprintf('============================================================\n\n');

fprintf('Tau PET vs RAW p-tau217\n');
fprintf(' Pearson r  = %.4f\n',rTauRaw);
fprintf(' Spearman rho = %.4f\n\n',rhoTauRaw);

fprintf('Tau PET vs LOG p-tau217\n');
fprintf(' Pearson r  = %.4f\n',rTauLog);
fprintf(' Spearman rho = %.4f\n\n',rhoTauLog);

fprintf('Absolute Pearson correlation:\n');
fprintf(' Raw p-tau217: %.4f\n',abs(rTauRaw));
fprintf(' Log p-tau217: %.4f\n',abs(rTauLog));
fprintf(' Difference (log - raw): %.4f\n',...
    abs(rTauLog)-abs(rTauRaw));

%% ------------------------------------------------------------
% Save outputs
% -------------------------------------------------------------

writetable(PearsonTable,...
    fullfile(outputDir,...
    'Step48A_Pearson_BiomarkerCorrelationMatrix.csv'),...
    'WriteRowNames',true);

writetable(SpearmanTable,...
    fullfile(outputDir,...
    'Step48A_Spearman_BiomarkerCorrelationMatrix.csv'),...
    'WriteRowNames',true);

writetable(PairwiseResults,...
    fullfile(outputDir,...
    'Step48A_Pairwise_BiomarkerCorrelations.csv'));

save(fullfile(outputDir,...
    'Step48A_BiomarkerCorrelations.mat'),...
    'X','biomarkerNames',...
    'Rpearson','Ppearson',...
    'Rspearman','Pspearman',...
    'PairwiseResults');

%% ------------------------------------------------------------
% Complete
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' STEP 48A COMPLETE\n');
fprintf('============================================================\n');

fprintf('Participant-level baseline analysis completed.\n');
fprintf('Expected participants: 143\n');
fprintf('Observed participants: %d\n',height(B));

fprintf('\nSaved outputs to:\n%s\n',outputDir);

fprintf('\nFiles created:\n');
fprintf('  Step48A_Pearson_BiomarkerCorrelationMatrix.csv\n');
fprintf('  Step48A_Spearman_BiomarkerCorrelationMatrix.csv\n');
fprintf('  Step48A_Pairwise_BiomarkerCorrelations.csv\n');
fprintf('  Step48A_BiomarkerCorrelations.mat\n');

fprintf('============================================================\n');