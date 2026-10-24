%% ================================================================
% Step38C_Final_Figure1_Reproducible.m
% ADNI REVISED COHORT — FIGURE 1
% VERSION 3: STRUCTURALLY CORRECTED LAYOUT
%
% Panel A: Cohort selection and exclusions
% Panel B: Study design and assessment timing
%
% Locked final cohort: N = 143
% CDR-SB observations: 379 (143 baseline + 236 post-index)
%
% Exports 600-dpi TIFF, vector PDF and PNG preview.
% Does not overwrite the approved 18 October reference figure.
%% ================================================================

clearvars;
clc;
close all;

%% DIRECTORIES



outputDir = 'D:\output';
figureDir = 'D:\output';

assert(isfolder(outputDir));
assert(isfolder(figureDir));

fprintf('\n============================================\n');
fprintf('STEP 38C — FIGURE 1 VERSION 3\n');
fprintf('============================================\n');

%% LOAD LOCKED DATA

F = readtable(fullfile(outputDir, ...
    'Step01_02_FTP_Cohort_Flow.csv'), ...
    'VariableNamingRule','preserve');

files = {
    'Step02_Earliest_Usable_FTP_Index_Cohort.csv'
    'Step03A_Amyloid_Matched_90Day_Cohort.csv'
    'Step03B_Amyloid_Positive_Cohort.csv'
    'Step04A_PTau217_90Day_Cohort.csv'
    'Step04B_Baseline_CDRSB_90Day_Cohort.csv'
    'Step05A_After_AdvancedDementia_Exclusion.csv'
    'Step05B_Final_143_Participant_Cohort.csv'
};

N = zeros(8,1);
N(1) = F.N_Reproduced(1);

for k = 1:7
    T = readtable(fullfile(outputDir,files{k}), ...
        'VariableNamingRule','preserve');
    N(k+1) = height(T);
end

T143 = readtable(fullfile(outputDir,files{7}), ...
    'VariableNamingRule','preserve');

T379 = readtable(fullfile(outputDir, ...
    'Step06_Final_Longitudinal_CDRSB_379_Observations.csv'), ...
    'VariableNamingRule','preserve');

expectedN = [1246;1237;1100;469;275;211;210;143];

assert(isequal(N,expectedN),'Cohort count mismatch.');

excluded = N(1:end-1)-N(2:end);

assert(isequal(excluded,[9;137;631;194;64;1;67]), ...
    'Exclusion count mismatch.');

assert(height(T379)==379);
assert(numel(unique(T379.RID))==143);

%% VERIFY TIMING

timingVariables = {
    'Amyloid_minus_FTP_days'
    'PTau217_minus_FTP_days'
    'CDRSB_minus_FTP_days'
};

medians = zeros(3,1);
quartiles = zeros(3,2);
ranges = zeros(3,2);

for k = 1:3

    x = double(T143.(timingVariables{k}));

    assert(numel(x)==143);
    assert(all(isfinite(x)));
    assert(all(abs(x)<=90));

    medians(k) = median(x);
    quartiles(k,:) = prctile(x,[25 75]);
    ranges(k,:) = [min(x) max(x)];

end

assert(isequal(medians,[-1;-2;-38]));
assert(isequal(ranges,[-35 84;-56 68;-89 5]));

fprintf('Cohort and timing verification: PASS\n');

%% FIGURE SETTINGS

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 13.5 12], ...
    'Renderer','painters');

set(fig,'DefaultAxesFontName','Arial');
set(fig,'DefaultTextFontName','Arial');

K = [0 0 0];
G = [0.35 0.35 0.35];
L = [0.96 0.96 0.96];
M = [0.90 0.90 0.90];

%% ================================================================
% PANEL A — COHORT SELECTION
%% ================================================================

axA = axes(fig,'Position',[0.02 0.34 0.96 0.64]);

axis(axA,[0 100 0 100]);
axis(axA,'off');
hold(axA,'on');

rectangle(axA,'Position',[0.5 0.5 99 99], ...
    'Curvature',0.01,'EdgeColor',K,'LineWidth',1.1);

text(axA,2,96,'A','FontSize',21,'FontWeight','bold');

% Main title occupies only the left/main column.
rectangle(axA,'Position',[14 91.5 50 7], ...
    'FaceColor',L,'EdgeColor',G,'Curvature',0.03);

text(axA,39,96, ...
    'ADNI COHORT SELECTION', ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',13);

text(axA,39,93.5, ...
    'For longitudinal CDR-SB analysis', ...
    'HorizontalAlignment','center','FontSize',10);

% Exclusion heading is in a separate column.
rectangle(axA,'Position',[68 91.5 28.5 7], ...
    'FaceColor',L,'EdgeColor',G,'Curvature',0.03);

text(axA,82.25,95, ...
    'Excluded at this step', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontWeight','bold','FontSize',10);

stageTitles = {
    'Participants with FTP tau PET'
    'Earliest usable FTP index scan'
    'Amyloid PET within ±90 days'
    'Amyloid-positive'
    'Fujirebio plasma p-tau217 within ±90 days'
    'Valid baseline CDR-SB within ±90 days'
    'Baseline CDR-SB < 9.5'
    'Final longitudinal cohort'
};

stageDetails = {
    'All available FTP scans'
    'Valid meta-temporal SUVR'
    'Any amyloid PET'
    'Berkeley amyloid status = 1'
    'Nearest valid plasma measurement'
    'Complete required baseline biomarkers'
    'Advanced baseline impairment excluded'
    '≥1 post-index CDR-SB, 3 months to 4 years'
};

exclusionDetails = {
    'No usable FTP index scan'
    'No amyloid PET within ±90 days'
    'Amyloid-negative or status unavailable'
    'No plasma p-tau217 within ±90 days'
    'No valid baseline CDR-SB within ±90 days'
    'Baseline CDR-SB ≥9.5'
    'No qualifying post-index CDR-SB'
};

yStage = [82 73.5 65 56.5 48 39.5 31 22.5];

mainX = 11;
mainW = 50;
mainH = 7;

exX = 68;
exW = 28.5;

for k = 1:8

    y = yStage(k);

    % Number box
    rectangle(axA,'Position',[5.5 y 3.8 mainH], ...
        'FaceColor',L,'EdgeColor',G,'Curvature',0.12);

    text(axA,7.4,y+3.5,num2str(k), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',12,'FontWeight','bold');

    % Cohort box
    rectangle(axA,'Position',[mainX y mainW mainH], ...
        'FaceColor',L,'EdgeColor',G,'Curvature',0.08);

    text(axA,mainX+1.2,y+4.8,stageTitles{k}, ...
        'FontSize',9.1,'FontWeight','bold', ...
        'Interpreter','none');

    text(axA,mainX+1.2,y+2.1,stageDetails{k}, ...
        'FontSize',8.1,'Interpreter','none');

    text(axA,mainX+mainW-1.3,y+4.8, ...
        sprintf('N = %s',commaNumber(N(k))), ...
        'HorizontalAlignment','right', ...
        'FontWeight','bold','FontSize',9.6);

    text(axA,mainX+mainW-1.3,y+2.1, ...
        'unique participants', ...
        'HorizontalAlignment','right','FontSize',7.8);

    % Cohort flow arrow
    if k<8

        xArrow = 36;

        plot(axA,[xArrow xArrow], ...
            [y yStage(k+1)+mainH], ...
            'Color',K,'LineWidth',0.9);

        plot(axA,xArrow,yStage(k+1)+mainH,'kv', ...
            'MarkerFaceColor',K,'MarkerSize',4);

    end

    % Exclusion box
    if k<=7

        rectangle(axA,'Position',[exX y exW mainH], ...
            'FaceColor',L,'EdgeColor',G,'Curvature',0.08);

        text(axA,exX+1.1,y+4.8, ...
            sprintf('Excluded: %s',commaNumber(excluded(k))), ...
            'FontSize',9.1,'FontWeight','bold');

        text(axA,exX+1.1,y+2.1,exclusionDetails{k}, ...
            'FontSize',8.1,'Interpreter','none');

        plot(axA,[mainX+mainW+0.8 exX-1.3], ...
            [y+mainH/2 y+mainH/2], ...
            '--','Color',G,'LineWidth',0.9);

        plot(axA,exX-1.3,y+mainH/2,'k>', ...
            'MarkerFaceColor',K,'MarkerSize',4);

    end

end

% Final dataset box
rectangle(axA,'Position',[13.5 4 48 12], ...
    'FaceColor',M,'EdgeColor',G,'Curvature',0.07);

text(axA,37.5,13.2,'Longitudinal dataset', ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',11);

text(axA,37.5,10.4,'N = 143 participants', ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',10);

text(axA,37.5,7.8,'379 CDR-SB observations', ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',10);

text(axA,37.5,5.6,'143 baseline + 236 post-index', ...
    'HorizontalAlignment','center','FontSize',9);

plot(axA,[36 36],[22.5 16],'k-','LineWidth',0.9);

%% ================================================================
% PANEL B — STUDY DESIGN AND TIMING
%
% Structural layout:
% Left: baseline assessment distributions
% Right: longitudinal follow-up
% Bottom: shared time axis
%% ================================================================

axB = axes(fig,'Position',[0.02 0.02 0.96 0.30]);

axis(axB,[0 100 0 100]);
axis(axB,'off');
hold(axB,'on');

rectangle(axB,'Position',[0.5 0.5 99 99], ...
    'Curvature',0.01,'EdgeColor',K,'LineWidth',1.1);

text(axB,2,94,'B','FontSize',21,'FontWeight','bold');

text(axB,50,94, ...
    'STUDY DESIGN AND TIMING OF ASSESSMENTS', ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',13);

text(axB,50,87, ...
    'Relative to index FTP tau PET', ...
    'HorizontalAlignment','center','FontSize',10);

%% Baseline timing region

% Baseline region occupies x = 3–46.
% Timing statistics occupy x = 27–46.
% Follow-up annotation occupies x = 50–96.

text(axB,22,79,'Index FTP tau PET (day 0)', ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',9);

rectangle(axB,'Position',[11 70 22 6], ...
    'FaceColor',M,'EdgeColor','none');

text(axB,22,73,'Baseline window: ±90 days', ...
    'HorizontalAlignment','center', ...
    'FontSize',8,'FontWeight','bold');

baselineRows = [60 48 36];

baselineLabels = {
    'Amyloid PET'
    'Plasma p-tau217'
    'Baseline CDR-SB'
};

% Local baseline-day scale, with day zero at x = 19.
baseX = @(d) 19 + d*(11/90);

plot(axB,[baseX(0) baseX(0)],[31 68], ...
    'k-','LineWidth',0.9);

for k = 1:3

    y = baselineRows(k);

    text(axB,2,y,baselineLabels{k}, ...
        'FontWeight','bold','FontSize',8.5, ...
        'VerticalAlignment','middle');

    % Full observed range
    plot(axB,baseX(ranges(k,:)),[y y], ...
        'k-','LineWidth',1);

    % Range caps
    for j = 1:2

        x = baseX(ranges(k,j));

        plot(axB,[x x],[y-1.3 y+1.3], ...
            'k-','LineWidth',0.8);

    end

    % Median
    symbols = {'o','s','d'};

    plot(axB,baseX(medians(k)),y, ...
        'Marker',symbols{k}, ...
        'MarkerFaceColor',K, ...
        'MarkerEdgeColor',K, ...
        'MarkerSize',5, ...
        'LineStyle','none');

    % Dedicated statistics region
    stats = sprintf( ...
        'Median %g d (IQR %g to %g)\nRange %g to %g d', ...
        medians(k),quartiles(k,1),quartiles(k,2), ...
        ranges(k,1),ranges(k,2));

    text(axB,33,y,stats, ...
        'FontSize',7.7, ...
        'VerticalAlignment','middle');

end

%% Longitudinal follow-up region

rectangle(axB,'Position',[51 37 45 31], ...
    'FaceColor',M,'EdgeColor','none');

text(axB,73.5,62, ...
    'Longitudinal CDR-SB assessments', ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',10);

text(axB,73.5,55, ...
    '3 months to 4 years after index FTP', ...
    'HorizontalAlignment','center','FontSize',9);

text(axB,73.5,48, ...
    '236 post-index observations', ...
    'HorizontalAlignment','center','FontSize',9);

text(axB,73.5,41, ...
    'Median last follow-up: 2.00 years', ...
    'HorizontalAlignment','center','FontSize',8.5);

%% Follow-up schematic arrow

plot(axB,[51 96],[28 28], ...
    'k-','LineWidth',1.2);

plot(axB,96,28,'k>', ...
    'MarkerFaceColor',K,'MarkerSize',5);

schematicX = [58 66 74 82 90];

for k = 1:numel(schematicX)

    plot(axB,schematicX(k),28,'ko', ...
        'MarkerFaceColor',K,'MarkerSize',4);

end

text(axB,73.5,22, ...
    'Schematic follow-up visits (not individual visit dates)', ...
    'HorizontalAlignment','center','FontSize',7.7);

%% Time-axis labels

text(axB,22,19,'Index FTP = day 0', ...
    'HorizontalAlignment','center','FontSize',8.5);

text(axB,51,16,'3 months', ...
    'HorizontalAlignment','center','FontSize',8.5);

text(axB,96,16,'4 years', ...
    'HorizontalAlignment','center','FontSize',8.5);

text(axB,50,6, ...
    'Assessment timing relative to index FTP tau PET', ...
    'HorizontalAlignment','center','FontSize',9.5);

%% EXPORT

tiffFile = fullfile(figureDir, ...
    'Figure1_Final_Step38C_600dpi.tif');

pdfFile = fullfile(figureDir, ...
    'Figure1_Final_Step38C_Vector.pdf');

pngFile = fullfile(figureDir, ...
    'Figure1_Final_Step38C_Preview.png');

exportgraphics(fig,tiffFile, ...
    'Resolution',600,'BackgroundColor','white');

exportgraphics(fig,pdfFile, ...
    'ContentType','vector','BackgroundColor','white');

exportgraphics(fig,pngFile, ...
    'Resolution',200,'BackgroundColor','white');

info = imfinfo(tiffFile);

fprintf('\nEXPORTED FIGURE\n');
fprintf('TIFF: %s\n',tiffFile);
fprintf('PDF: %s\n',pdfFile);
fprintf('PNG: %s\n',pngFile);

fprintf('TIFF resolution: %.0f x %.0f dpi\n', ...
    info.XResolution,info.YResolution);

fprintf('\nSTEP 38C VERSION 3 COMPLETE\n');
fprintf('Data verification: PASS\n');
fprintf('Visual QC: PENDING\n');
fprintf('Approved reference: UNCHANGED\n');

%% LOCAL FUNCTION

function s = commaNumber(x)

    s = sprintf('%.0f',x);

    if x<0
        s = ['-' commaNumber(abs(x))];
        return;
    end

    for k = length(s)-2:-3:2
        s = [s(1:k-1) ',' s(k:end)];
    end

end