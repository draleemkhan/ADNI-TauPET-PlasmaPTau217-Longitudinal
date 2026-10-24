%% STEP 59 - DEFINITIVE APOE MATCHING TO LOCKED N=143 (MATLAB R2025a)
% Requires Step58 QC PASS. Does not overwrite any Step01-Step58 outputs.
% Primary source: APOERES; C2N fallback only where APOERES unavailable.
% Conflict detection is repeated independently from Step58 audit records.
clear; clc;

outDir = 'D:\output';
cohortFile = fullfile(outDir,'Step05B_Final_143_Participant_Cohort.csv');
sourceFile = fullfile(outDir,'Step58_APOE_AllSourceRecords.csv');
auditFile = fullfile(outDir,'Step58_APOE_CohortRID_Audit.csv');

for f = [string(cohortFile),string(sourceFile),string(auditFile)]
    assert(isfile(f),'Required input missing: %s',f);
end
logFile = fullfile(outDir,'Step59_APOE_Final143_Matching_Log.txt');
if isfile(logFile), error('Step59 output exists; archive it before rerunning: %s',logFile); end
diary(logFile); diary on;
try
    fprintf('STEP 59 | DEFINITIVE APOE MATCHING\n');
    C = readtable(cohortFile,'TextType','string','VariableNamingRule','preserve');
    S = readtable(sourceFile,'TextType','string','VariableNamingRule','preserve');
    A = readtable(auditFile,'TextType','string','VariableNamingRule','preserve');

    assert(all(ismember({'RID'},C.Properties.VariableNames)),'Cohort RID missing');
    assert(all(ismember({'RID','Source','Genotype','ValidGenotype','ValidRID'},S.Properties.VariableNames)), ...
        'Step58 source schema mismatch');
    assert(all(ismember({'RID','Conflict'},A.Properties.VariableNames)),'Step58 audit schema mismatch');
    rid = str2double(string(C.RID));
    assert(height(C)==143 && all(isfinite(rid)) && numel(unique(rid))==143, ...
        'Locked N=143 cohort failed validation');
    assert(height(A)==143 && ~any(logical(A.Conflict)), ...
        'Step58 audit not full-coverage or contains conflicts');
    S.RID = str2double(string(S.RID));
    S.Source = string(S.Source);
    S.Genotype = string(S.Genotype);
    valid = ismember(S.RID,rid) & logical(S.ValidRID) & logical(S.ValidGenotype);
    Q = S(valid,:);
    n=height(C);
    G = strings(n,1); G(:)=missing;
    Source = strings(n,1); Source(:)=missing;
    N_APOERES=zeros(n,1); N_C2N=zeros(n,1);
    for k=1:n
        rows = Q(Q.RID==rid(k),:);
        assert(~isempty(rows),'Missing genotype for RID %g',rid(k));
        ug=unique(rows.Genotype);
        assert(numel(ug)==1,'Discordant genotypes for RID %g',rid(k));
        ap=rows(rows.Source=="APOERES",:);
        c2=rows(rows.Source=="C2N",:);
        N_APOERES(k)=height(ap); N_C2N(k)=height(c2);
        if ~isempty(ap)
            G(k)=ap.Genotype(1);
            Source(k)="APOERES";
        else
            assert(~isempty(c2),'No C2N fallback for RID %g',rid(k));
            G(k)=c2.Genotype(1);
            Source(k)="C2N fallback";
        end
    end
    allowed=["2/2","2/3","2/4","3/3","3/4","4/4"];
    assert(all(ismember(G,allowed)),'Invalid selected genotype');
    alleleCount=double(contains(G,"4"));
    alleleCount(G=="4/4")=2;
    carrier=double(alleleCount>=1);
    % Preserve locked cohort columns unchanged; append new uniquely named fields.
    Final=C;
    assert(~any(ismember({'APOE59_Genotype','APOE59_Source','APOE59_AlleleCount','APOE59_Carrier'}, ...
        Final.Properties.VariableNames)),'Step59 variable collision');
    Final.APOE59_Genotype=G;
    Final.APOE59_Source=Source;
    Final.APOE59_AlleleCount=alleleCount;
    Final.APOE59_Carrier=carrier;
    Final.APOE59_N_APOERES=N_APOERES;
    Final.APOE59_N_C2N=N_C2N;
    assert(height(Final)==143 && numel(unique(rid))==143,'Cohort changed');
    assert(all(~ismissing(G)),'Missing genotype after merge');
    Summary=table( ...
        ["Locked cohort N";"APOE coverage";"APOERES selected";"C2N fallback selected"; ...
         "APOE4 carriers";"APOE4 noncarriers"], ...
        [n;sum(~ismissing(G));sum(Source=="APOERES");sum(Source=="C2N fallback"); ...
         sum(carrier==1);sum(carrier==0)], ...
        'VariableNames',{'Metric','Value'});
    [gg,gn]=findgroups(G);
    Distribution=table(gn,splitapply(@numel,G,gg),'VariableNames',{'Genotype','N'});
    disp(Summary); disp(Distribution);
    % Refuse overwrite of any Step59 output.
    f1=fullfile(outDir,'Step59_Final143_With_APOE.csv');
    f2=fullfile(outDir,'Step59_APOE_QC_Summary.csv');
    f3=fullfile(outDir,'Step59_APOE_Genotype_Distribution.csv');
    for f=[string(f1),string(f2),string(f3)]
        assert(~isfile(f),'Output already exists; refusing overwrite: %s',f);
    end
    writetable(Final,f1);
    writetable(Summary,f2);
    writetable(Distribution,f3);
    fprintf('STEP 59 QC PASS. Saved outputs without modifying earlier steps.\n');
catch ME
    fprintf(2,'STEP 59 FAILED: %s\n',getReport(ME,'extended','hyperlinks','off'));
    diary off; rethrow(ME);
end
diary off;
