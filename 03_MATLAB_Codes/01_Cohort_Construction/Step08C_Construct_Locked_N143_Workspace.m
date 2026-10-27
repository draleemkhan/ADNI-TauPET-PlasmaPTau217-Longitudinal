%% Step08C_Construct_Locked_N143_Workspace.m
% ADNI revised cohort - locked longitudinal workspace reconstruction
% Input: Step08_ModelReady_379_Observations.csv
% Output: ADNI_REVISED_COHORT_LOCKED_N143.mat
% The original locked workspace was created interactively; this script reconstructs it from Step08 output.
% Verified identical to the original on 26 Oct 2026 (all 21 variables; exact MATLAB table equality).
% Output contains participant-level ADNI data and must remain local.

clc;

root = 'D:\';
inputFile = fullfile(root,'output','Step08_ModelReady_379_Observations.csv');
outputDir = fullfile(root,'output');
outputFile = fullfile(outputDir,'ADNI_REVISED_COHORT_LOCKED_N143.mat');

assert(isfile(inputFile),'Step08 input CSV not found.');
assert(~isfile(outputFile),'Locked workspace already exists. Refusing to overwrite.');

%% Read upstream model-ready data
T = readtable(inputFile,'VariableNamingRule','preserve');

%% Verify upstream structure
assert(height(T)==379,'Unexpected observation count.');
assert(numel(unique(T.RID))==143,'Unexpected participant count.');
assert(width(T)==19,'Unexpected Step08 variable count.');

%% Construct categorical modelling variables
T.RIDcat = categorical(T.RID);
T.SexModel = categorical(string(T.Sex),["Male","Female"]);

%% Verify reconstructed structure
assert(height(T)==379,'Reconstructed row count mismatch.');
assert(width(T)==21,'Reconstructed variable count mismatch.');
assert(~any(isundefined(T.RIDcat)),'Undefined participant category.');
assert(~any(isundefined(T.SexModel)),'Undefined sex category.');

%% Preserve expected workspace variable name
LockedLongitudinal = T;

%% Save locked workspace
save(outputFile,'LockedLongitudinal','-v7.3');

%% Verify saved file can be reloaded
S = load(outputFile,'LockedLongitudinal');
assert(isequaln(S.LockedLongitudinal,LockedLongitudinal), ...
    'Locked workspace reload verification failed.');

fprintf('\n============================================\n');
fprintf(' STEP08C - LOCKED WORKSPACE CONSTRUCTION\n');
fprintf('============================================\n');
fprintf('Participants: %d\n',numel(unique(LockedLongitudinal.RID)));
fprintf('Observations: %d\n',height(LockedLongitudinal));
fprintf('Variables: %d\n',width(LockedLongitudinal));
fprintf('Locked workspace reload verification: PASS\n');
fprintf('Locked workspace saved to:\n%s\n',outputFile);
fprintf('============================================\n');
