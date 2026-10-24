%% STEP 60: FREEZE FULL-COHORT APOE LONGITUDINAL DATA (MATLAB R2025a)
% Continues locked ADNI pipeline. Does not fit models or overwrite old files.
clear; clc;

outDir = 'D:\output';
cohortFile = fullfile(outDir,'Step05B_Final_143_Participant_Cohort.csv');
genoFile = fullfile(outDir,'Step59_Final143_With_APOE.csv');
modelFile = fullfile(outDir,'Step08_ModelReady_379_Observations.csv');
logFile = fullfile(outDir,'Step60_APOE_Longitudinal_Freeze_Log.txt');
outputFile = fullfile(outDir,'Step60_ModelReady_379_With_APOE.csv');
participantFile = fullfile(outDir,'Step60_APOE_Participant_Audit.csv');
summaryFile = fullfile(outDir,'Step60_APOE_QC_Summary.csv');
for f = [string(cohortFile),string(genoFile),string(modelFile)]
    assert(isfile(f),'Required input not found: %s',f);
end
for f = [string(logFile),string(outputFile),string(participantFile),string(summaryFile)]
    assert(~isfile(f),'Refusing to overwrite existing output: %s',f);
end
diary(logFile); diary on;
try
    fprintf('STEP 60 | APOE LONGITUDINAL DATA FREEZE\n');
    C=readtable(cohortFile,'TextType','string','VariableNamingRule','preserve');
    A=readtable(genoFile,'TextType','string','VariableNamingRule','preserve');
    D=readtable(modelFile,'TextType','string','VariableNamingRule','preserve');
    requiredA={'RID','APOE59_Genotype','APOE59_Source','APOE59_AlleleCount','APOE59_Carrier'};
    requiredD={'RID','CDRSB','TimeYears','Tau_z','Centiloid_z','PTau217_z','Age','Sex','Education'};
    assert(all(ismember(requiredA,A.Properties.VariableNames)),'Step59 schema mismatch');
    assert(all(ismember(requiredD,D.Properties.VariableNames)),'Step08 schema mismatch');
    assert(ismember('RID',C.Properties.VariableNames),'Step05B missing RID');
    ridC=str2double(string(C.RID)); ridA=str2double(string(A.RID)); ridD=str2double(string(D.RID));
    assert(height(C)==143 && numel(unique(ridC))==143 && all(isfinite(ridC)), 'Locked cohort failed QC');
    assert(height(A)==143 && numel(unique(ridA))==143 && isequal(sort(ridA),sort(ridC)), 'Step59 RIDs differ from locked cohort');
    assert(height(D)==379 && numel(unique(ridD))==143 && isequal(sort(unique(ridD)),sort(ridC)), 'Step08 observations/RIDs differ');
    assert(~any(ismember(requiredA(2:end),D.Properties.VariableNames)),'APOE field name collision');
    [tf,loc]=ismember(ridD,ridA);
    assert(all(tf),'Unmatched Step08 RID');
    M=D;
    for j=2:numel(requiredA)
        nm=requiredA{j};
        M.(nm)=A.(nm)(loc,:);
    end
    assert(height(M)==379 && isequaln(M(:,D.Properties.VariableNames),D), ...
        'Original model-ready rows or values were altered');
    assert(all(ismember(string(M.APOE59_Genotype),["2/2","2/3","2/4","3/3","3/4","4/4"])), 'Invalid genotype');
    assert(all(ismember(M.APOE59_Carrier,[0,1])),'Invalid carrier status');
    assert(all(M.APOE59_Carrier==double(M.APOE59_AlleleCount>=1)), 'Carrier/count disagreement');
    assert(all(ismember(M.APOE59_AlleleCount,[0,1,2])),'Invalid allele count');
    P=A(:,requiredA);
    assert(sum(P.APOE59_Carrier==1)==95 && sum(P.APOE59_Carrier==0)==48, 'Carrier counts differ from Step59');
    assert(sum(string(P.APOE59_Source)=="APOERES")==142 && ...
           sum(string(P.APOE59_Source)=="C2N fallback")==1,'Provenance differs from Step59');
    [~,locP]=ismember(ridC,ridA);
    P=P(locP,:); % participant audit in locked cohort order
    obsCount=zeros(height(P),1);
    for k=1:height(P), obsCount(k)=sum(ridD==ridC(k)); end
    P.N_LongitudinalObservations=obsCount;
    assert(all(obsCount>=2) && sum(obsCount)==379,'Longitudinal counts inconsistent');
    Metric=["Participants";"Observations";"Unique RIDs";"APOE4 carriers";"APOE4 noncarriers"; ...
        "APOERES selected";"C2N fallback selected";"Participants with >=2 observations"];
    Value=[height(P);height(M);numel(unique(ridD));sum(P.APOE59_Carrier==1); ...
        sum(P.APOE59_Carrier==0);sum(string(P.APOE59_Source)=="APOERES"); ...
        sum(string(P.APOE59_Source)=="C2N fallback");sum(obsCount>=2)];
    QC=table(Metric,Value);
    disp(QC);
    writetable(M,outputFile);
    writetable(P,participantFile);
    writetable(QC,summaryFile);
    fprintf('STEP 60 QC PASS. Original Step08 rows/values preserved; genotype appended.\n');
    fprintf('Saved: %s\n',outputFile);
catch ME
    fprintf(2,'STEP 60 FAILED: %s\n',getReport(ME,'extended','hyperlinks','off'));
    diary off; rethrow(ME);
end
diary off;
