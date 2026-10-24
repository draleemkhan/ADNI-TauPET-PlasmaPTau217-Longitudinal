%% STEP 62: FULL-COHORT APOE4-ADJUSTED REML LME (MATLAB R2025a)
% Continuation of locked ADNI pipeline, following Step61 QC PASS.
% Fits SAME N=143 / 379 rows in both models; no restandardization.
% A: exact Step61 primary formula. B: A + APOE4 carrier + carrier:time.
% Does not overwrite any earlier outputs.
clear; clc;

outDir = 'D:\output';
inFile = fullfile(outDir,'Step60_ModelReady_379_With_APOE.csv');
step61File = fullfile(outDir,'Step61_PrimaryModel_Reproduction.csv');
originalFile = fullfile(outDir,'Step08_ModelReady_379_Observations.csv');
files = {fullfile(outDir,'Step62_APOE_Adjusted_LME_Results.csv'), ...
    fullfile(outDir,'Step62_APOE_Adjusted_LME_Comparison.csv'), ...
    fullfile(outDir,'Step62_APOE_Adjusted_LME_Coefficients.csv'), ...
    fullfile(outDir,'Step62_APOE_Adjusted_LME.mat'), ...
    fullfile(outDir,'Step62_APOE_Adjusted_LME_Log.txt')};
assert(isfile(inFile) && isfile(step61File) && isfile(originalFile), ...
    'Missing Step60, Step61, or Step08 input.');
for k=1:numel(files)
    assert(~isfile(files{k}),'Refusing to overwrite: %s',files{k});
end
D = readtable(inFile,'VariableNamingRule','preserve');
O = readtable(originalFile,'VariableNamingRule','preserve');
R61 = readtable(step61File,'VariableNamingRule','preserve');
required = {'RID','CDRSB','TimeYears','Tau_z','Centiloid_z', ...
    'PTau217_z','Age','Sex','Education','APOE59_Carrier','APOE59_AlleleCount'};
assert(all(ismember(required,D.Properties.VariableNames)), ...
    'Step60 is missing required model/APOE columns.');
assert(height(D)==379 && height(O)==379 && height(R61)==1, ...
    'Unexpected longitudinal or Step61 row count.');
assert(numel(unique(D.RID))==143 && isequaln(D.RID,O.RID), ...
    'RID cohort size/order differs from locked Step08.');
for k=1:9
    v=required{k};
    assert(ismember(v,O.Properties.VariableNames),'Step08 missing %s',v);
    if isnumeric(D.(v)) && isnumeric(O.(v))
        assert(isequaln(D.(v),O.(v)),'Original field changed: %s',v);
    else
        assert(isequaln(string(D.(v)),string(O.(v))), ...
            'Original field changed: %s',v);
    end
end
for k=2:numel(required)
    x=D.(required{k});
    if isnumeric(x)
        assert(all(isfinite(x)),'Nonfinite values in %s',required{k});
    else
        assert(~any(ismissing(x)),'Missing values in %s',required{k});
    end
end
assert(all(ismember(D.APOE59_Carrier,[0,1])) && ...
    all(ismember(D.APOE59_AlleleCount,[0,1,2])) && ...
    all(D.APOE59_Carrier==double(D.APOE59_AlleleCount>0)), ...
    'APOE coding inconsistent.');
[~,ia]=unique(D.RID,'stable');
assert(sum(D.APOE59_Carrier(ia)==1)==95 && ...
    sum(D.APOE59_Carrier(ia)==0)==48,'APOE participant distribution mismatch.');
for i=1:numel(ia)
    idx=D.RID==D.RID(ia(i));
    assert(numel(unique(D.APOE59_Carrier(idx)))==1, ...
        'Within-RID APOE inconsistency.');
end
D.RIDcat=categorical(D.RID);
% Step56B / Step61 primary formula reproduced verbatim.
formulaA = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
    'Age + Sex + Education + ' ...
    'Tau_z:TimeYears + Centiloid_z:TimeYears + ' ...
    'PTau217_z:TimeYears + ' ...
    '(1|RIDcat) + (-1 + TimeYears|RIDcat)'];
% Step56B adjusted formula, replacing legacy APOE4_Carrier with verified APOE59_Carrier.
formulaB = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
    'Age + Sex + Education + APOE59_Carrier + ' ...
    'Tau_z:TimeYears + Centiloid_z:TimeYears + ' ...
    'PTau217_z:TimeYears + APOE59_Carrier:TimeYears + ' ...
    '(1|RIDcat) + (-1 + TimeYears|RIDcat)'];
fprintf('STEP 62 | FULL-COHORT APOE4-ADJUSTED REML LME\n');
mdlA=fitlme(D,formulaA,'FitMethod','REML');
CA=mdlA.Coefficients;
namesA=string(CA.Name);
rowA=find(contains(namesA,'Tau_z') & contains(namesA,'TimeYears'));
assert(numel(rowA)==1,'Reference tau x time coefficient ambiguous.');
ciA=coefCI(mdlA);
betaA=CA.Estimate(rowA); seA=CA.SE(rowA);
loA=ciA(rowA,1); hiA=ciA(rowA,2); pA=CA.pValue(rowA);
assert(abs(betaA-R61.TauTime_Beta(1))<1e-6 && ...
    abs(seA-R61.TauTime_SE(1))<1e-6 && ...
    abs(pA-R61.TauTime_PValue(1))<1e-6, ...
    'Reference model does not reproduce Step61. STOP.');
mdlB=fitlme(D,formulaB,'FitMethod','REML');
CB=mdlB.Coefficients;
namesB=string(CB.Name);
rowB=find(contains(namesB,'Tau_z') & contains(namesB,'TimeYears'));
rowAPOE=find(contains(namesB,'APOE59_Carrier') & contains(namesB,'TimeYears'));
assert(numel(rowB)==1 && numel(rowAPOE)==1, ...
    'Adjusted tau or APOE interaction coefficient ambiguous.');
ciB=coefCI(mdlB);
betaB=CB.Estimate(rowB); seB=CB.SE(rowB);
loB=ciB(rowB,1); hiB=ciB(rowB,2); pB=CB.pValue(rowB);
apoeBeta=CB.Estimate(rowAPOE); apoeSE=CB.SE(rowAPOE);
apoeLo=ciB(rowAPOE,1); apoeHi=ciB(rowAPOE,2);
apoeP=CB.pValue(rowAPOE);
Model=["Primary reference (no APOE)";"APOE4-adjusted"];
NParticipants=[143;143]; NObservations=[379;379];
TauTime_Beta=[betaA;betaB]; TauTime_SE=[seA;seB];
TauTime_CI_Lower=[loA;loB]; TauTime_CI_Upper=[hiA;hiB];
TauTime_PValue=[pA;pB];
Results=table(Model,NParticipants,NObservations,TauTime_Beta, ...
    TauTime_SE,TauTime_CI_Lower,TauTime_CI_Upper,TauTime_PValue);
DeltaBeta=betaB-betaA;
PercentChange=100*DeltaBeta/betaA;
Comparison=table(betaA,betaB,DeltaBeta,PercentChange, ...
    apoeBeta,apoeSE,apoeLo,apoeHi,apoeP, ...
    'VariableNames',{'ReferenceTauTimeBeta','AdjustedTauTimeBeta', ...
    'DeltaBeta','PercentChange','APOE4_Time_Beta','APOE4_Time_SE', ...
    'APOE4_Time_CI_Lower','APOE4_Time_CI_Upper','APOE4_Time_PValue'});
% Handle fitlme coefficient outputs returned as table or legacy dataset.
if isa(CA,'dataset'), CA=dataset2table(CA); end
if isa(CB,'dataset'), CB=dataset2table(CB); end
assert(istable(CA) && istable(CB),'Unexpected coefficient object type.');
CA.Model=repmat("Primary reference",height(CA),1);
CB.Model=repmat("APOE4-adjusted",height(CB),1);
Coefficients=[CA;CB];
fprintf('\n'); disp(Results); disp(Comparison);
% Save only after both models, checks and export table preparation succeed.
writetable(Results,files{1});
writetable(Comparison,files{2});
writetable(Coefficients,files{3});
save(files{4},'mdlA','mdlB','Results','Comparison','Coefficients', ...
    'formulaA','formulaB','inFile','step61File');
fid=fopen(files{5},'w');
assert(fid>0,'Unable to write Step62 log.');
fprintf(fid,'STEP 62: FULL-COHORT APOE4-ADJUSTED REML LME\n');
fprintf(fid,'Input: %s\nStep61 reference: %s\n',inFile,step61File);
fprintf(fid,'N=143, observations=379, carriers=95, noncarriers=48\n');
fprintf(fid,'Reference formula: %s\nAdjusted formula: %s\n',formulaA,formulaB);
fprintf(fid,'Reference tau*time: beta=%.9g SE=%.9g CI=[%.9g,%.9g] p=%.9g\n',betaA,seA,loA,hiA,pA);
fprintf(fid,'Adjusted tau*time: beta=%.9g SE=%.9g CI=[%.9g,%.9g] p=%.9g\n',betaB,seB,loB,hiB,pB);
fprintf(fid,'Delta beta=%.9g; percent change=%.9g%%\n',DeltaBeta,PercentChange);
fprintf(fid,'APOE4*time: beta=%.9g SE=%.9g CI=[%.9g,%.9g] p=%.9g\n',apoeBeta,apoeSE,apoeLo,apoeHi,apoeP);
fprintf(fid,'REML used for both models. No REML likelihood-ratio comparison performed.\n');
fclose(fid);
fprintf('STEP 62 QC PASS. Reference reproduces Step61; adjusted model fitted.\n');
fprintf('Saved: %s\n',files{1});
