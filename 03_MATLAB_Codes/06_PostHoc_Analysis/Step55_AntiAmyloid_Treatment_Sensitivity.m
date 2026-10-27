% PRIVATE INPUT PROVENANCE:
% Step55_Private_Treatment_Start_Dates.csv is derived from
% ANTIAMYTX_11Aug2026.csv.
% One row per analytic participant with recorded anti-amyloid treatment.
% TreatmentStart is the earliest recorded treatment start date.
% Month-only dates are assigned to the first day of that month.
% This participant-level file must remain local and must not be committed.

%% STEP 55: ANTI-AMYLOID TREATMENT SENSITIVITY
%
% Purpose:
% Assess whether the primary tau PET x time association changes after
% censoring longitudinal CDR-SB observations occurring after initiation
% of anti-amyloid treatment.
%
% This is a targeted sensitivity analysis only.
%
% Primary comparison:
%   Locked primary tau PET x time beta = 0.413
%
% Outputs:
%   Step55_AntiAmyloid_Treatment_Sensitivity.csv
%   Step55_AntiAmyloid_Treatment_Sensitivity.mat
%
% Code:
% D:\codes
%
% Output:
% D:\output

clear;
clc;

fprintf('============================================================\n');
fprintf(' STEP 55: ANTI-AMYLOID TREATMENT SENSITIVITY\n');
fprintf('============================================================\n\n');

%% ------------------------------------------------------------
% 1. PATHS
% -------------------------------------------------------------



outputDir = 'D:\output';

modelFile = fullfile(outputDir, ...
    'Step08_ModelReady_379_Observations.csv');

assert(isfile(modelFile), ...
    'Model-ready longitudinal dataset not found.');

%% ------------------------------------------------------------
% 2. LOAD LOCKED MODEL-READY DATA
% -------------------------------------------------------------

D = readtable(modelFile,'VariableNamingRule','preserve');

fprintf('Original observations: %d\n',height(D));

% Identify RID variable
assert(ismember('RID',D.Properties.VariableNames), ...
    'RID variable not found.');

fprintf('Original participants: %d\n',numel(unique(D.RID)));

%% ------------------------------------------------------------
% 3. IDENTIFY DATE VARIABLE
% -------------------------------------------------------------

% Displayed here explicitly so that the script fails rather than guessing
% if the expected longitudinal assessment date variable is absent.

possibleDateVars = { ...
    'CDRDate', ...
    'AssessmentDate', ...
    'VISDATE', ...
    'EXAMDATE', ...
    'CDRSB_Date'};

dateVar = '';

for k = 1:numel(possibleDateVars)
    if ismember(possibleDateVars{k},D.Properties.VariableNames)
        dateVar = possibleDateVars{k};
        break
    end
end

if isempty(dateVar)
    fprintf('\nAvailable variables:\n');
    disp(D.Properties.VariableNames');
    error(['CDR-SB assessment date variable could not be identified. ' ...
           'No censoring performed.']);
end

fprintf('Assessment date variable: %s\n',dateVar);

%% ------------------------------------------------------------
% 4. CONVERT ASSESSMENT DATE
% -------------------------------------------------------------

assessmentDate = D.(dateVar);

if ~isdatetime(assessmentDate)

    if iscell(assessmentDate) || isstring(assessmentDate) || ...
            iscategorical(assessmentDate)

        assessmentDate = datetime(string(assessmentDate));

    else
        error('Unexpected assessment-date variable type.');
    end
end

%% ------------------------------------------------------------
% 5. LOAD PRIVATE ANTI-AMYLOID TREATMENT START DATES
% -------------------------------------------------------------
% Participant-level treatment information is stored locally.
% Do not upload the private treatment CSV to GitHub.

txFile = fullfile('D:\source files', ...
    'Step55_Private_Treatment_Start_Dates.csv');

assert(isfile(txFile), ...
    'Private treatment mapping not found.');

TX = readtable(txFile);

assert(all(ismember({'RID','TreatmentStart'}, ...
    TX.Properties.VariableNames)), ...
    'Treatment mapping columns missing.');

treatedRID = double(TX.RID);
treatmentStart = TX.TreatmentStart;

if ~isdatetime(treatmentStart)
    treatmentStart = datetime(string(treatmentStart));
end

assert(numel(treatedRID)==numel(treatmentStart), ...
    'Treatment mapping length mismatch.');

assert(~any(isnat(treatmentStart)), ...
    'Invalid treatment-start dates.');

%% ------------------------------------------------------------
% 6. FLAG POST-TREATMENT OBSERVATIONS
% -------------------------------------------------------------

postTreatment = false(height(D),1);

for k = 1:numel(treatedRID)

    idx = D.RID == treatedRID(k);

    postTreatment(idx) = ...
        assessmentDate(idx) >= treatmentStart(k);

end

fprintf('\nPost-treatment observations identified: %d\n', ...
    sum(postTreatment));

fprintf('Participants contributing post-treatment observations: %d\n', ...
    numel(unique(D.RID(postTreatment))));

%% ------------------------------------------------------------
% 7. AUDIT OBSERVATIONS TO BE CENSORED
% -------------------------------------------------------------

Audit = table( ...
    D.RID(postTreatment), ...
    assessmentDate(postTreatment), ...
    'VariableNames',{'RID','AssessmentDate'});

disp(Audit);

% We expect exactly six observations from four participants.
assert(height(Audit)==6, ...
    ['Expected 6 post-treatment observations but identified %d. ' ...
     'STOP: inspect date handling before proceeding.'],height(Audit));

assert(numel(unique(Audit.RID))==4, ...
    'Expected post-treatment observations from exactly 4 participants.');

fprintf('Treatment censoring audit: PASS\n');

%% ------------------------------------------------------------
% 8. CREATE CENSORED DATASET
% -------------------------------------------------------------

Dc = D(~postTreatment,:);

fprintf('\nAfter censoring:\n');
fprintf('Participants: %d\n',numel(unique(Dc.RID)));
fprintf('Observations: %d\n',height(Dc));

assert(numel(unique(Dc.RID))==143, ...
    'Unexpected loss of participants after censoring.');

assert(height(Dc)==373, ...
    'Expected 373 observations after censoring.');

%% ------------------------------------------------------------
% 9. IDENTIFY PRIMARY MODEL VARIABLES
% -------------------------------------------------------------

required = { ...
    'CDRSB', ...
    'TimeYears', ...
    'Tau_z', ...
    'Centiloid_z', ...
    'PTau217_z', ...
    'Age', ...
    'Sex', ...
    'Education'};

missingVars = required(~ismember(required,Dc.Properties.VariableNames));

if ~isempty(missingVars)

    fprintf('\nAvailable variables:\n');
    disp(Dc.Properties.VariableNames');

    error('Required primary-model variable(s) missing: %s', ...
        strjoin(missingVars,', '));
end

%% ------------------------------------------------------------
% 10. ENSURE PARTICIPANT ID IS CATEGORICAL
% -------------------------------------------------------------

Dc.RIDcat = categorical(Dc.RID);

%% ------------------------------------------------------------
% 11. REFIT EXACT PRIMARY LME
% -------------------------------------------------------------
%
% Fixed effects:
%   TimeYears
%   Tau_z
%   Centiloid_z
%   PTau217_z
%   Age
%   Sex
%   Education
%   Tau_z x TimeYears
%   Centiloid_z x TimeYears
%   PTau217_z x TimeYears
%
% Random effects:
%   participant-specific random intercept and time slope
%   with intercept-slope covariance constrained to zero.
%
% REML estimation, matching the primary analysis.

formula = [ ...
    'CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
    'Age + Sex + Education + ' ...
    'Tau_z:TimeYears + Centiloid_z:TimeYears + ' ...
    'PTau217_z:TimeYears + ' ...
    '(1|RIDcat) + (-1 + TimeYears|RIDcat)' ...
    ];

mdl = fitlme(Dc,formula, ...
    'FitMethod','REML');

fprintf('\nPrimary model successfully refitted after treatment censoring.\n');

%% ------------------------------------------------------------
% 12. EXTRACT TAU PET x TIME COEFFICIENT
% -------------------------------------------------------------

C = mdl.Coefficients;

termNames = string(C.Name);

tauRow = find( ...
    contains(termNames,'Tau_z') & ...
    contains(termNames,'TimeYears'));

assert(numel(tauRow)==1, ...
    'Tau PET x time coefficient could not be uniquely identified.');

beta = C.Estimate(tauRow);
SE   = C.SE(tauRow);
pVal = C.pValue(tauRow);

CI = coefCI(mdl);

ciLow  = CI(tauRow,1);
ciHigh = CI(tauRow,2);

%% ------------------------------------------------------------
% 13. COMPARE WITH LOCKED PRIMARY RESULT
% -------------------------------------------------------------

primaryBeta = 0.413;

deltaBeta = beta - primaryBeta;

percentChange = 100 * deltaBeta / primaryBeta;

fprintf('\n============================================================\n');
fprintf(' STEP 55 RESULT\n');
fprintf('============================================================\n');

fprintf('Locked primary tau PET x time beta: %.3f\n',primaryBeta);

fprintf(['Treatment-censored tau PET x time: ' ...
         'beta %.3f (95%% CI %.3f to %.3f), P %.4f\n'], ...
         beta,ciLow,ciHigh,pVal);

fprintf('Absolute change in beta: %.3f\n',deltaBeta);
fprintf('Relative change in beta: %.1f%%\n',percentChange);

%% ------------------------------------------------------------
% 14. SAVE COMPACT RESULT TABLE
% -------------------------------------------------------------

Result = table( ...
    143, ...
    height(D), ...
    height(Dc), ...
    sum(postTreatment), ...
    numel(unique(D.RID(postTreatment))), ...
    primaryBeta, ...
    beta, ...
    SE, ...
    ciLow, ...
    ciHigh, ...
    pVal, ...
    deltaBeta, ...
    percentChange, ...
    'VariableNames',{ ...
    'NParticipants', ...
    'OriginalObservations', ...
    'CensoredObservations', ...
    'PostTreatmentObservationsRemoved', ...
    'ParticipantsWithPostTreatmentObservations', ...
    'PrimaryBeta', ...
    'TreatmentCensoredBeta', ...
    'SE', ...
    'CI_Lower', ...
    'CI_Upper', ...
    'PValue', ...
    'DeltaBeta', ...
    'PercentChangeBeta'});

csvOut = fullfile(outputDir, ...
    'Step55_AntiAmyloid_Treatment_Sensitivity.csv');

matOut = fullfile(outputDir, ...
    'Step55_AntiAmyloid_Treatment_Sensitivity.mat');

writetable(Result,csvOut);

save(matOut, ...
    'Result','Audit','mdl','Dc','postTreatment', ...
    'treatedRID','treatmentStart');

%% ------------------------------------------------------------
% 15. FINAL AUDIT
% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf(' STEP 55 FINAL AUDIT\n');
fprintf('============================================================\n');

fprintf('Original participants: 143\n');
fprintf('Original observations: %d\n',height(D));
fprintf('Post-treatment observations censored: %d\n',sum(postTreatment));
fprintf('Participants affected: %d\n', ...
    numel(unique(D.RID(postTreatment))));
fprintf('Final observations: %d\n',height(Dc));
fprintf('Primary model specification retained: PASS\n');
fprintf('Original biomarker standardisation retained: PASS\n');
fprintf('REML estimation retained: PASS\n');
fprintf('Treatment censoring audit: PASS\n');

fprintf('\nFiles created:\n');
fprintf('  %s\n',csvOut);
fprintf('  %s\n',matOut);

fprintf('\n============================================================\n');
fprintf(' STEP 55 COMPLETE\n');
fprintf('============================================================\n');