%% Step08_Demographics_and_ModelReady_Dataset.m
% ADNI Tau PET revised cohort reconstruction
%
% PURPOSE:
% 1) Add Age, Sex, and Education to the final 143-participant cohort.
% 2) Use PTDEMOG as the authoritative demographic source.
% 3) Construct a 379-row model-ready longitudinal dataset by joining
%    participant-level baseline variables to Step06 repeated CDR-SB rows.
%
% Demographic rules:
% - Sex: PTGENDER (1=Male, 2=Female)
% - Education: PTEDUCAT (years)
% - Birth year: PTDOBYY
% - Age at index FTP = elapsed years from PTDOBYY date to index FTP date.
%   PTDOBYY is stored by ADNI as a year-level date (e.g. 1931-01-01).
% - When multiple demographic records exist, use the record closest in
%   absolute time to index FTP among records with valid sex/education/DOB.
%
% This script DOES NOT fit the mixed-effects model yet.

clear; clc; close all;

%% Locked paths
sourceDir = 'D:\source files';
outputDir = 'D:\output';

participantFile = fullfile(outputDir,'Step07_Final143_With_Baseline_Diagnosis.csv');
longFile        = fullfile(outputDir,'Step06_Final_Longitudinal_CDRSB_379_Observations.csv');
demogFile       = fullfile(sourceDir,'PTDEMOG_11Aug2026.csv');

assert(isfile(participantFile),'Missing Step07 participant file.');
assert(isfile(longFile),'Missing Step06 longitudinal file.');
assert(isfile(demogFile),'Missing PTDEMOG source file.');

C = readtable(participantFile,'VariableNamingRule','preserve');
L = readtable(longFile,'VariableNamingRule','preserve');
D = readtable(demogFile,'VariableNamingRule','preserve');

fprintf('\n====================================================\n');
fprintf(' STEP 8 - DEMOGRAPHICS + MODEL-READY DATASET\n');
fprintf('====================================================\n');
fprintf('Participant cohort: %d\n',height(C));
fprintf('Repeated CDR-SB rows: %d\n',height(L));

%% Required PTDEMOG fields
req = {'RID','VISDATE','PTGENDER','PTDOBYY','PTEDUCAT'};
missing = setdiff(req,D.Properties.VariableNames);
if ~isempty(missing)
    error('PTDEMOG missing required variables: %s',strjoin(missing,', '));
end

%% Index dates
if ismember('IndexFTP_Date',C.Properties.VariableNames)
    idxDate = parseDateFlexible(C.IndexFTP_Date);
elseif ismember('SCANDATE',C.Properties.VariableNames)
    idxDate = parseDateFlexible(C.SCANDATE);
else
    error('Index FTP date not found in participant cohort.');
end

%% Demographic source conversion
dRID = double(D.RID);
dVisit = parseDateFlexible(D.VISDATE);
dDOB = parseDateFlexible(D.PTDOBYY);

if isnumeric(D.PTGENDER)
    dSex = double(D.PTGENDER);
else
    dSex = str2double(string(D.PTGENDER));
end

if isnumeric(D.PTEDUCAT)
    dEduc = double(D.PTEDUCAT);
else
    dEduc = str2double(string(D.PTEDUCAT));
end

validSex  = ismember(dSex,[1 2]);
validEduc = isfinite(dEduc) & dEduc >= 0 & dEduc <= 30;
validDOB  = ~isnat(dDOB);

%% Match one demographic record per participant
n = height(C);

SexCode = nan(n,1);
Education = nan(n,1);
DOB = NaT(n,1);
DemogVisitDate = NaT(n,1);
DemogGapDays = nan(n,1);

for i = 1:n
    rid = double(C.RID(i));
    cand = find(dRID==rid & validSex & validEduc & validDOB);

    if isempty(cand)
        continue
    end

    % Prefer closest dated demographic record to index FTP.
    dated = cand(~isnat(dVisit(cand)));

    if ~isempty(dated)
        gaps = days(dVisit(dated)-idxDate(i));
        [~,j] = min(abs(gaps));
        k = dated(j);
        DemogGapDays(i) = gaps(j);
        DemogVisitDate(i) = dVisit(k);
    else
        k = cand(1);
    end

    SexCode(i) = dSex(k);
    Education(i) = dEduc(k);
    DOB(i) = dDOB(k);
end

%% Age at index FTP
Age = years(idxDate - DOB);

Sex = strings(n,1);
Sex(SexCode==1) = "Male";
Sex(SexCode==2) = "Female";
Sex(~ismember(SexCode,[1 2])) = "Missing";

%% Add to participant cohort
C.Age = Age;
C.SexCode = SexCode;
C.Sex = Sex;
C.Education = Education;
C.DOB_YearDate = DOB;
C.Demog_VisitDate = DemogVisitDate;
C.Demog_minus_FTP_days = DemogGapDays;

%% Audit completeness
fprintf('\nDEMOGRAPHIC COMPLETENESS\n');
fprintf('Age available:       %d / %d\n',sum(isfinite(C.Age)),n);
fprintf('Sex available:       %d / %d\n',sum(ismember(C.Sex,["Male","Female"])),n);
fprintf('Education available: %d / %d\n',sum(isfinite(C.Education)),n);

if all(isfinite(C.Age)) && all(ismember(C.Sex,["Male","Female"])) && all(isfinite(C.Education))
    fprintf('DEMOGRAPHIC COMPLETENESS: PASS\n');
else
    fprintf('DEMOGRAPHIC COMPLETENESS: REVIEW\n');
end

%% Descriptive audit
fprintf('\nDEMOGRAPHIC SUMMARY\n');
fprintf('Age: mean %.2f | SD %.2f | range %.2f-%.2f years\n', ...
    mean(C.Age,'omitnan'),std(C.Age,'omitnan'),min(C.Age),max(C.Age));
fprintf('Education: mean %.2f | SD %.2f | range %.0f-%.0f years\n', ...
    mean(C.Education,'omitnan'),std(C.Education,'omitnan'), ...
    min(C.Education),max(C.Education));
fprintf('Male: %d | Female: %d\n',sum(C.Sex=="Male"),sum(C.Sex=="Female"));

%% Verify required baseline biomarker fields exist
needed = {'RID','META_TEMPORAL_SUVR','Centiloid','PTau217_Fujirebio', ...
          'Baseline_Diagnosis','Age','Sex','Education'};
missing2 = setdiff(needed,C.Properties.VariableNames);
if ~isempty(missing2)
    error('Participant cohort missing model variables: %s',strjoin(missing2,', '));
end

%% Build participant covariate table
P = C(:,needed);
P.Properties.VariableNames{'META_TEMPORAL_SUVR'} = 'TauMetaSUVR';

% One row per RID check
if height(P) ~= numel(unique(P.RID))
    error('Participant covariate table is not one row per RID.');
end

%% Join with 379 repeated CDR-SB rows
ModelData = innerjoin(L,P,'Keys','RID');

fprintf('\nMODEL-READY LONGITUDINAL DATASET\n');
fprintf('Rows after join: %d [expected 379]\n',height(ModelData));
fprintf('Unique participants: %d [expected 143]\n',numel(unique(ModelData.RID)));

if height(ModelData)==379 && numel(unique(ModelData.RID))==143
    fprintf('MODEL-READY JOIN: PASS\n');
else
    fprintf('MODEL-READY JOIN: REVIEW\n');
end

%% Convert diagnosis and sex to categorical with explicit references
ModelData.Diagnosis = categorical(string(ModelData.Baseline_Diagnosis), ...
    ["Cognitively normal","MCI","Dementia"]);
ModelData.Sex_cat = categorical(string(ModelData.Sex),["Male","Female"]);

%% Z-standardize continuous baseline predictors
% Candidate model specified by investigator:
% Tau, Centiloid, p-tau217 are z-standardized.
ModelData.Tau_z = zscore(double(ModelData.TauMetaSUVR));
ModelData.Centiloid_z = zscore(double(ModelData.Centiloid));
ModelData.PTau217_z = zscore(double(ModelData.PTau217_Fujirebio));

% Save standardization parameters for exact reproducibility
StdVars = ["TauMetaSUVR";"Centiloid";"PTau217_Fujirebio"];
Means = [mean(ModelData.TauMetaSUVR,'omitnan'); ...
         mean(ModelData.Centiloid,'omitnan'); ...
         mean(ModelData.PTau217_Fujirebio,'omitnan')];
SDs = [std(ModelData.TauMetaSUVR,'omitnan'); ...
       std(ModelData.Centiloid,'omitnan'); ...
       std(ModelData.PTau217_Fujirebio,'omitnan')];

% Because baseline values repeat across outcome rows, compute participant-level
% parameters separately as well; these are the parameters to use for the model.
ParticipantMeans = [mean(P.TauMetaSUVR,'omitnan'); ...
                    mean(P.Centiloid,'omitnan'); ...
                    mean(P.PTau217_Fujirebio,'omitnan')];
ParticipantSDs = [std(P.TauMetaSUVR,'omitnan'); ...
                  std(P.Centiloid,'omitnan'); ...
                  std(P.PTau217_Fujirebio,'omitnan')];

% IMPORTANT correction: z-standardization must use participant-level baseline
% distributions, not repeated-row weighted distributions.
ModelData.Tau_z = (double(ModelData.TauMetaSUVR)-ParticipantMeans(1))/ParticipantSDs(1);
ModelData.Centiloid_z = (double(ModelData.Centiloid)-ParticipantMeans(2))/ParticipantSDs(2);
ModelData.PTau217_z = (double(ModelData.PTau217_Fujirebio)-ParticipantMeans(3))/ParticipantSDs(3);

Standardization = table(StdVars,ParticipantMeans,ParticipantSDs, ...
    'VariableNames',{'Variable','ParticipantMean','ParticipantSD'});

fprintf('\nPARTICIPANT-LEVEL STANDARDIZATION PARAMETERS\n');
disp(Standardization);

%% Save outputs
writetable(C,fullfile(outputDir,'Step08_Final143_With_Demographics.csv'));
writetable(ModelData,fullfile(outputDir,'Step08_ModelReady_379_Observations.csv'));
writetable(Standardization,fullfile(outputDir,'Step08_Standardization_Parameters.csv'));

fprintf('\nOutputs saved to:\n%s\n',outputDir);
fprintf('====================================================\n');

%% Local date parser
function d = parseDateFlexible(x)
    if isdatetime(x)
        d=x; d.Format='yyyy-MM-dd'; return
    end

    s=strtrim(string(x));
    d=NaT(size(s));

    fmts={'yyyy-MM-dd','yyyy-MM-dd HH:mm:ss','MM/dd/yyyy','MM/yyyy', ...
          'dd-MMM-yyyy','dd/MM/yyyy'};

    rem=~ismissing(s) & strlength(s)>0;

    for k=1:numel(fmts)
        if ~any(rem), break; end
        try
            tmp=datetime(s(rem),'InputFormat',fmts{k});
            ok=~isnat(tmp);
            idx=find(rem);
            d(idx(ok))=tmp(ok);
            rem(idx(ok))=false;
        catch
        end
    end

    if any(rem)
        try
            tmp=datetime(s(rem));
            idx=find(rem);
            ok=~isnat(tmp);
            d(idx(ok))=tmp(ok);
        catch
        end
    end

    d.Format='yyyy-MM-dd';
end
