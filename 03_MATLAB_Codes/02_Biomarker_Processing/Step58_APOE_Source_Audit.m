%% STEP 58 - APOE SOURCE AUDIT (MATLAB R2025a)
% Continuation of locked Steps 56A/56B and 57A-57C.
% Reads source CSVs and locked cohort; DOES NOT change any legacy outputs.
% Output is source-level QC only; Step 59 will perform the definitive merge.
clear; clc;

inDir = 'D:\source files';
outDir = 'D:\output';
if ~isfolder(outDir), mkdir(outDir); end
apoeFile = fullfile(inDir,'APOERES_18Oct2026.csv');
c2nFile = fullfile(inDir,'C2N_PRECIVITYAD2_PLASMA_11Aug2026.csv');
cohortFile = fullfile(outDir,'Step05B_Final_143_Participant_Cohort.csv');
assert(isfile(apoeFile),'Missing APOERES: %s',apoeFile);
assert(isfile(c2nFile),'Missing C2N: %s',c2nFile);
assert(isfile(cohortFile),'Missing locked cohort: %s',cohortFile);
logFile = fullfile(outDir,'Step58_APOE_Source_Audit_Log.txt');
if isfile(logFile), delete(logFile); end
diary(logFile); diary on;
try
    fprintf('STEP 58 | APOE SOURCE AUDIT\n');
    fprintf('APOERES: %s\nC2N: %s\nCohort: %s\n',apoeFile,c2nFile,cohortFile);
    A = readtable(apoeFile,'TextType','string','VariableNamingRule','preserve');
    B = readtable(c2nFile,'TextType','string','VariableNamingRule','preserve');
    C = readtable(cohortFile,'VariableNamingRule','preserve');
    assert(all(ismember({'RID','GENOTYPE'},A.Properties.VariableNames)),'APOERES schema mismatch');
    assert(all(ismember({'RID','APOE_C2N'},B.Properties.VariableNames)),'C2N schema mismatch');
    assert(ismember('RID',C.Properties.VariableNames),'Cohort RID missing');
    rC = str2double(string(C.RID));
    assert(height(C)==143 && all(isfinite(rC)) && numel(unique(rC))==143,'Locked cohort QC failed');
    rA = str2double(string(A.RID)); rB = str2double(string(B.RID));
    gA = standardizeAPOE(string(A.GENOTYPE));
    gB = standardizeAPOE(string(B.APOE_C2N));
    SA = table(rA,repmat("APOERES",height(A),1),string(A.GENOTYPE),gA,...
        'VariableNames',{'RID','Source','RawGenotype','Genotype'});
    SB = table(rB,repmat("C2N",height(B),1),string(B.APOE_C2N),gB,...
        'VariableNames',{'RID','Source','RawGenotype','Genotype'});
    S = [SA;SB];
    S.InLockedCohort = ismember(S.RID,rC);
    S.ValidRID = isfinite(S.RID) & S.RID>0 & S.RID==floor(S.RID);
    S.ValidGenotype = ~ismissing(S.Genotype);
    S.InvalidNonmissing = ~S.ValidGenotype & ~ismissing(strtrim(S.RawGenotype)) & strlength(strtrim(S.RawGenotype))>0 & ...
        ~ismember(upper(strtrim(S.RawGenotype)),["NA","N/A","NAN","NULL",".","<MISSING>"]);
    % Audit every within-source and between-source record; do not silently resolve conflicts.
    Q = S(S.InLockedCohort & S.ValidRID & S.ValidGenotype,:);
    ridList = unique(Q.RID);
    perRID = table('Size',[numel(ridList),6], ...
        'VariableTypes',{'double','double','double','string','logical','logical'}, ...
        'VariableNames',{'RID','N_APOERES','N_C2N','Genotypes','Conflict','HasAPOERES'});
    for k=1:numel(ridList)
        rid=ridList(k); rows=Q(Q.RID==rid,:); ug=unique(rows.Genotype);
        perRID.RID(k)=rid;
        perRID.N_APOERES(k)=sum(rows.Source=="APOERES");
        perRID.N_C2N(k)=sum(rows.Source=="C2N");
        perRID.Genotypes(k)=strjoin(ug,' | ');
        perRID.Conflict(k)=numel(ug)>1;
        perRID.HasAPOERES(k)=any(rows.Source=="APOERES");
    end
    missingRID = setdiff(rC,ridList);
    conflicts = perRID(perRID.Conflict,:);
    invalid = S(S.InvalidNonmissing | (~S.ValidRID),:);
    % Source-level and cohort-specific genotype counts, no unverified expected counts hard-coded.
    fprintf('\nRaw records: APOERES=%d, C2N=%d\n',height(A),height(B));
    fprintf('Locked cohort N=%d\n',numel(rC));
    fprintf('Cohort RIDs with >=1 valid genotype=%d\n',height(perRID));
    fprintf('Cohort RIDs with no valid genotype=%d\n',numel(missingRID));
    fprintf('Cohort RIDs with discordant genotypes=%d\n',height(conflicts));
    fprintf('Cohort RIDs represented by APOERES=%d\n',sum(perRID.HasAPOERES));
    fprintf('Cohort RIDs represented by C2N only=%d\n',sum(~perRID.HasAPOERES));
    fprintf('Invalid nonmissing or invalid RID source records=%d\n',height(invalid));
    fprintf('Records with valid genotype in cohort=%d\n',height(Q));
    writetable(S,fullfile(outDir,'Step58_APOE_AllSourceRecords.csv'));
    writetable(perRID,fullfile(outDir,'Step58_APOE_CohortRID_Audit.csv'));
    writetable(conflicts,fullfile(outDir,'Step58_APOE_Conflicts.csv'));
    writetable(invalid,fullfile(outDir,'Step58_APOE_InvalidRecords.csv'));
    writetable(table(missingRID,'VariableNames',{'RID'}),fullfile(outDir,'Step58_APOE_MissingRIDs.csv'));
    if ~isempty(conflicts)
        fprintf('QC HOLD: conflicting genotypes require review before Step 59.\n');
    elseif isempty(missingRID)
        fprintf('QC PASS: full cohort coverage and no genotype conflicts.\n');
    else
        fprintf('QC PARTIAL: missing cohort genotype records require review.\n');
    end
    fprintf('STEP 58 COMPLETED. No prior outputs changed.\n');
catch ME
    fprintf(2,'STEP 58 FAILED: %s\n',getReport(ME,'extended','hyperlinks','off'));
    diary off; rethrow(ME);
end
diary off;

function out = standardizeAPOE(raw)
    out = strings(size(raw)); out(:)=missing;
    for j=1:numel(raw)
        s=upper(strtrim(raw(j)));
        if ismissing(s) || ismember(s,["","NA","N/A","NAN","NULL",".","<MISSING>"])
            continue;
        end
        s=replace(s,["EPSILON","Ε","Ɛ","E"," ","-","_"],["","","","","","",""]);
        s=replace(s,["\\","|",",",";"],["/","/","/","/"]);
        if ~contains(s,'/') && strlength(s)==2 && all(ismember(char(s),'234'))
            s=extractBetween(s,1,1)+"/"+extractBetween(s,2,2);
        end
        parts=split(s,'/');
        if numel(parts)~=2, continue; end
        a=str2double(parts);
        if any(~ismember(a,[2,3,4])), continue; end
        a=sort(a);
        out(j)=string(a(1))+"/"+string(a(2));
    end
end
