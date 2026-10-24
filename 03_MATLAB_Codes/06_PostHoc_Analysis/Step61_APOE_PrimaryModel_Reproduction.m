%% STEP 61: REPRODUCE LOCKED PRIMARY LME BEFORE APOE ADJUSTMENT
% Reads frozen Step60 data; DOES NOT include APOE in the model.
% Preserves original standardized biomarkers and all 379 observations.
% MATLAB R2025a / Statistics and Machine Learning Toolbox.
clear; clc;

outDir = 'D:\output';
inputFile = fullfile(outDir,'Step60_ModelReady_379_With_APOE.csv');
originalFile = fullfile(outDir,'Step08_ModelReady_379_Observations.csv');
assert(isfile(inputFile),'Missing Step60 frozen data: %s',inputFile);
assert(isfile(originalFile),'Missing original Step08 data: %s',originalFile);
outputs = {fullfile(outDir,'Step61_PrimaryModel_Reproduction.csv'), ...
           fullfile(outDir,'Step61_PrimaryModel_Coefficients.csv'), ...
           fullfile(outDir,'Step61_PrimaryModel_Reproduction.mat'), ...
           fullfile(outDir,'Step61_PrimaryModel_Reproduction_Log.txt')};
for k=1:numel(outputs)
    assert(~isfile(outputs{k}),'Refusing to overwrite existing output: %s',outputs{k});
end
D = readtable(inputFile,'VariableNamingRule','preserve');
O = readtable(originalFile,'VariableNamingRule','preserve');
required = {'RID','CDRSB','TimeYears','Tau_z','Centiloid_z','PTau217_z','Age','Sex','Education'};
assert(all(ismember(required,D.Properties.VariableNames)), ...
    'Missing required primary-model variables.');
assert(height(D)==379 && height(O)==379,'Expected 379 longitudinal observations.');
assert(numel(unique(D.RID))==143 && numel(unique(O.RID))==143, ...
    'Expected 143 unique participants.');
assert(isequaln(D.RID,O.RID),'Step60 RID order differs from Step08.');
% Confirm exact source values for the variables used by the model.
for k=1:numel(required)
    v=required{k};
    assert(ismember(v,O.Properties.VariableNames),'Step08 missing %s',v);
    x=D.(v); y=O.(v);
    if isnumeric(x) && isnumeric(y)
        assert(isequaln(x,y),'Step60 changed original variable %s',v);
    else
        assert(isequaln(string(x),string(y)), ...
            'Step60 changed original variable %s',v);
    end
end
% Explicit complete-case verification: no silently dropped model rows.
for k=2:numel(required)
    v=required{k}; x=D.(v);
    if isnumeric(x)
        assert(all(isfinite(x)),'Missing/nonfinite values in %s',v);
    else
        assert(~any(ismissing(x)),'Missing values in %s',v);
    end
end
D.RIDcat = categorical(D.RID);
% Formula reproduced verbatim from Step56B, Model A.
formula = ['CDRSB ~ TimeYears + Tau_z + Centiloid_z + PTau217_z + ' ...
    'Age + Sex + Education + ' ...
    'Tau_z:TimeYears + Centiloid_z:TimeYears + ' ...
    'PTau217_z:TimeYears + ' ...
    '(1|RIDcat) + (-1 + TimeYears|RIDcat)'];
mdl = fitlme(D,formula,'FitMethod','REML');
coef = mdl.Coefficients;
names = string(coef.Name);
tauRow = find(contains(names,'Tau_z') & contains(names,'TimeYears'));
assert(numel(tauRow)==1,'Tau PET x time term not uniquely identified.');
ci = coefCI(mdl);
beta = coef.Estimate(tauRow);
se = coef.SE(tauRow);
lo = ci(tauRow,1);
hi = ci(tauRow,2);
p = coef.pValue(tauRow);
% Historical rounded primary result: beta 0.413, 95%% CI 0.087-0.739, p .013.
% Compare within rounding precision, not exact binary floating-point equality.
betaMatches = abs(beta-0.413)<0.001;
ciMatches = abs(lo-0.087)<0.002 && abs(hi-0.739)<0.002;
pMatches = abs(p-0.013)<0.002;
Result = table(143,379,beta,se,lo,hi,p,betaMatches,ciMatches,pMatches, ...
    'VariableNames',{'NParticipants','NObservations','TauTime_Beta', ...
    'TauTime_SE','TauTime_CI_Lower','TauTime_CI_Upper','TauTime_PValue', ...
    'BetaMatchesPublished','CIMatchesPublished','PMatchesPublished'});
logID=fopen(outputs{4},'w');
assert(logID>0,'Could not create Step61 log.');
fprintf('STEP 61 | PRIMARY LME REPRODUCTION (NO APOE TERMS)\n');
disp(Result);
fprintf(logID,'Input: %s\nOriginal: %s\nFormula: %s\n',inputFile,originalFile,formula);
fprintf(logID,'N=143, observations=379\nBeta=%.9f, SE=%.9f, CI=[%.9f, %.9f], p=%.9g\n',beta,se,lo,hi,p);
fprintf(logID,'Matches historical rounded results: beta=%d, CI=%d, p=%d\n', ...
    betaMatches,ciMatches,pMatches);
fclose(logID);
% Fail safely if historical result is not reproduced: do not label PASS.
if ~(betaMatches && ciMatches && pMatches)
    warning('STEP61:ReproductionMismatch', ...
        'Primary result differs from historical rounded values. Inspect model and log before Step62.');
    fprintf('STEP 61 QC REVIEW REQUIRED. No result files written. Log saved.\n');
else
    writetable(Result,outputs{1});
    writetable(dataset2table(coef),outputs{2});
    save(outputs{3},'mdl','Result','formula','required','inputFile','originalFile');
    fprintf('STEP 61 QC PASS. Primary LME reproduced without APOE adjustment.\n');
    fprintf('Saved: %s\n',outputs{1});
end
