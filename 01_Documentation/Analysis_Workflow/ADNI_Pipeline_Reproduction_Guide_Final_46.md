# ADNI longitudinal tau PET study — pipeline reproduction guide

This guide describes how an authorised ADNI researcher can reproduce the study's cohort construction, longitudinal models, and sensitivity analyses using the accompanying MATLAB and R scripts. The repository provides **analysis code and input-data documentation**, not participant-level ADNI data. It is a methodological guide, **not a claim that a clean-environment end-to-end rerun has been independently completed**.

## 1. Obtain the inputs

Obtain the required ADNI releases through the applicable ADNI data-access process. Consult `01_Documentation/Data_Access/ADNI_Source_Data_Inventory_14.csv` for the exact filenames and purposes. The principal analytical inputs are UC Berkeley **non-partial-volume-corrected** 6-mm flortaucipir tau PET, UC Berkeley 6-mm amyloid PET, Fujirebio plasma p-tau217, CDR, DXSUM, PTDEMOG, APOERES, C2N PrecivityAD2 plasma and anti-amyloid treatment records. The five PET metadata files in the 14-file inventory are documented as potential supporting inputs; verify whether individual scripts require them. Do not commit restricted ADNI files or participant-level derived tables to a public repository.

**Source clarification:** The source-lineage audit independently reproduced all 1,237 index participants and their meta-temporal SUVRs from `UCBERKELEY_TAU_6MM_11Aug2026.csv`, the UC Berkeley non-partial-volume-corrected (NPVC) tau PET dataset. This verified NPVC dataset was used for cohort construction and baseline tau PET quantification. The partial-volume-corrected dataset (`UCBERKELEY_TAUPVC_6MM`) was not used.

## 2. Configure the environment

- MATLAB **R2025a**, with the toolboxes/functions used by the scripts (including `fitlme`), and R with **geepack** for GEE sensitivity analyses.
- **Agreed Windows reproduction directories:**
  - `D:\source files` — the 14 authorised ADNI source CSV files listed in `01_Documentation/Data_Access/ADNI_Source_Data_Inventory_14.csv`.
  - `D:\codes` — locally downloaded MATLAB analysis scripts.
  - `D:\output` — writable destination for intermediate CSVs, model results, tables and figures (TIFF/PDF/PNG). Figures are saved directly here, **without a separate figures subfolder**.
- **GitHub portability corrections:** The repository contains **48 MATLAB scripts** and **2 R scripts**. The identified historical absolute-directory references in GitHub copies were corrected to use `D:\source files` for authorised ADNI source CSVs and `D:\output` for intermediate files and generated results. The corresponding GitHub SHA-256 hashes were refreshed; original research scripts and their source hashes were not changed. These corrections do not establish that every script has been successfully executed in a clean environment.
- **Repository layout versus runtime layout:** The GitHub repository organises MATLAB scripts under `03_MATLAB_Codes/` in functional subfolders and stores the 14-file inventory under `01_Documentation/Data_Access/`. These repository paths are **not** the same as the agreed local `D:\codes` directory. Do not assume that flattening the scripts into `D:\codes` preserves relative-path lookups.
- **Known path dependency requiring confirmation:** `Step11B_Broad_ADNI_Registry_Activity_Audit.m` locates `ADNI_Source_Data_Inventory_14.csv` using a repository-relative path. Before executing it from `D:\codes`, either retain the expected repository-relative directory structure or verify an explicitly approved local inventory-path adjustment. **The flattened `D:\codes` setup has not yet been confirmed to satisfy this lookup.** Do not alter the SHA-256-locked script merely to make the guide appear turnkey.
- Source datasets must remain access-controlled. Do not publish `D:\source files` or participant-level outputs from `D:\output`.
- Preserve source filenames, variables, dates and units, and keep participant identifier `RID` consistent across tables. Review each script's actual path settings and input/output dependencies before execution; the collection is not yet a unified, parameterised one-command pipeline. Execute in dependency order.

## 3. Reconstruct the analytical cohort

1. **Index tau PET:** Select the earliest [18F]flortaucipir scan per RID with a valid meta-temporal SUVR from the verified non-PVC Berkeley dataset. Expected index cohort: **1,237**.
2. **Amyloid PET:** Match the nearest Berkeley amyloid PET within **±90 days** of index tau PET, then retain `AMYLOID_STATUS == 1`. Expected counts: **1,100** matched, **469** amyloid-positive. Use the corrected `Step03_Amyloid_Match_and_Positivity.m`; the historical benchmark of 1,101 is obsolete.

**Amyloid PET tracer eligibility:** Initial screening permitted [¹⁸F]florbetapir, [¹⁸F]florbetaben, and [¹⁸F]NAV4694. Following application of all prespecified eligibility criteria, the final analytical cohort of 143 participants contained 74 florbetapir and 69 florbetaben examinations, with no NAV4694 examinations. NAV4694 was not retrospectively excluded from the initial screening criteria.
3. **Plasma and baseline clinical status:** Select the valid Fujirebio `pT217_F` measurement nearest index tau PET within ±90 days and the nearest valid baseline CDR-SB assessment within ±90 days. Exclude baseline CDR-SB **≥9.5**.
4. **Longitudinal qualification:** Require at least one post-index CDR-SB assessment between **3 months and 4 years** after index tau PET. Retain eligible visits individually; do not average repeated observations. Date each observation relative to the **index tau PET**, including baseline CDR-SB at its actual signed assessment time.
5. **Diagnosis and demographics:** Assign diagnosis nearest index tau PET within ±90 days (CN/MCI/dementia). Obtain age, sex and education from demographic records, choosing the valid record nearest index tau PET where multiple records exist.

**Expected sequential cohort counts:** `1,246 → 1,237 → 1,100 → 469 → 275 → 211 → 210 → 143`. Final cohort: **143 participants** (33 CN, 75 MCI, 35 dementia), contributing **379 CDR-SB observations**. Check these against the script-generated flow and final model-ready table; they are reproduction checkpoints, not filters to hard-code.

**Step08C — Locked workspace reconstruction:** Run `Step08C_Construct_Locked_N143_Workspace.m` immediately after Step08 and before any downstream scripts requiring `ADNI_REVISED_COHORT_LOCKED_N143.mat`. Step08C reconstructs the locked N = 143 analytical workspace from the Step08 model-ready output.

## 4. Prepare biomarkers and time

Use baseline meta-temporal tau PET SUVR, Berkeley Centiloids and Fujirebio p-tau217 as continuous variables. Compute **z scores using one baseline value per participant in the final 143-person cohort**, and retain those original scaling constants in sensitivity subsets and bootstrap analyses. Reported SDs: tau PET **0.393 SUVR**, amyloid **35.88 Centiloids**, p-tau217 **0.357 pg/mL**. Express follow-up time continuously in years from index tau PET. Model sex categorically (male reference).

## 5. Fit the primary longitudinal model

Outcome: repeated CDR-SB. Fixed effects: time, standardised tau PET, Centiloids and p-tau217, age, sex, education, and **each biomarker × time** interaction. Include participant-specific random intercepts and random linear-time slopes with **zero intercept–slope covariance**. Fit the primary model using **restricted maximum likelihood (REML)**. The prespecified main parameter is **tau PET × time**, interpreted as the difference in annual CDR-SB change per 1-SD higher baseline tau PET, conditional on the other covariates. Do **not** include baseline CDR-SB as a primary adjustment covariate; it is part of the repeated outcome.

## 6. Compare incremental model information

Fit nested mixed models on **identical observations**, with the same random-effects structure, in this order: (a) Centiloids + demographics + Centiloid × time; (b) add p-tau217 and p-tau217 × time; (c) add tau PET and tau PET × time. Use **maximum likelihood (ML)**, not REML, for fixed-effect nested-model likelihood-ratio tests and AIC/BIC comparisons. Each biomarker block adds **two fixed-effect terms**. Calculate marginal and conditional R² using the study's variance decomposition for diagonal random effects. Assess the tau PET increment in marginal R² with **2,000 participant-cluster bootstrap resamples**, retaining full participant trajectories and original z-score scaling. These measures reflect **in-sample fit**, not externally validated prediction.

## 7. Reproduce robustness analyses

Run the relevant scripts for residual/influence diagnostics; alternative random-effects structures; quadratic time and tau × time²; unweighted and inverse-probability-weighted GEE (independence/exchangeable correlation, robust standard errors); diagnosis-adjusted and stage-specific models; **2,000-cluster** bootstrap for tau × time; timing restrictions (±60 and ±30 days); baseline clinical-severity adjustment; biomarker correlations and fixed-effect collinearity; longitudinal-qualification/attrition assessment, inverse-probability weighting and weight truncation; and the assumption-based MNAR tipping-point sensitivity. Maintain the same underlying model specification and scaling where the Methods require it.

Reproduce separately labelled **post hoc** analyses: reverse biomarker entry order; diagnosis × biomarker × time interactions; whole-cohort and MCI **1-df slope-specific** likelihood-ratio tests (distinct from 2-df biomarker-block tests); log-transformed p-tau217; anti-amyloid treatment censoring; and APOE ε4 adjustment. APOE uses APOERES preferentially with C2N fallback. Do not treat exploratory P values as multiplicity-adjusted confirmatory tests.

### Attrition-audit intermediate-file dependencies (2026-10-24 update)

The final **47-script GitHub MATLAB collection includes** the following three source-integrity-verified scripts, added earlier to document and complete the intermediate-file chain for longitudinal-qualification and attrition audits (they are **already counted** in the 44 scripts):

1. `Step11_Dropout_Censoring_Audit.m` — audits dropout and censoring among the baseline-eligible participants, using the upstream cohort-construction output, including `Step05A_After_AdvancedDementia_Exclusion.csv`.
2. `Step11B_Broad_ADNI_Registry_Activity_Audit.m` — audits subsequent activity across the documented 14-file ADNI source inventory (`01_Documentation/Data_Access/ADNI_Source_Data_Inventory_14.csv`) and generates `Step11B_BroadRegistry_ParticipantAudit.csv`.
3. `Step11C_Retained_vs_NonRetained_Attrition_Comparison_v3.m` — reads `Step11B_BroadRegistry_ParticipantAudit.csv` and demographic/clinical source files to compare retained and nonretained participants, including a subgroup with at least four years of potential follow-up opportunity.

**Execution guidance:** Generate the required upstream cohort outputs first. Run the Step 11 and Step 11B audits before Step 11C; **Step 11B's participant-audit CSV is a confirmed direct input to Step 11C**. This recommended order does not imply that Step 11B directly consumes a Step 11 output. Check each script's remaining input/output dependencies and paths in the final dependency-chain audit. Intermediate participant-level CSVs are local research outputs and must not be published in the repository.

## 8. Validate outputs and report transparently

Verify the cohort-flow counts, final **143 RIDs/379 observations**, biomarker standardisation, primary model coefficients, nested-model degrees of freedom, and figures/tables against the manuscript. The previously reported primary tau PET × time coefficient is approximately **0.413** (95% CI **0.087–0.739**, P≈**.013**); use this as a checkpoint, **not** a fixed target. Preserve generated outputs and any discrepancies. Distinguish independently rerun results from comparisons against saved summaries. `Source_SHA256.csv` lists 50 scripts (48 MATLAB and two R), recording the SHA-256 hash of each original source script, where applicable, and its GitHub copy. Step50B (post hoc slope-equality test) was added during manuscript revision and has no earlier source version. All GitHub-copy hashes matched the local files before commit (Windows `certutil`). Hash verification establishes file identity, not successful execution or reproduction of analytical results.

## 9. Scope and data-sharing limitations

All participant-level ADNI data must be acquired independently under the applicable ADNI conditions. The repository can share code, variable definitions, workflow instructions and non-identifying summary outputs. The source-code manifest establishes copy integrity, **not** proof of complete end-to-end execution on another system. Check each script for hard-coded paths, intermediate-file dependencies and local software assumptions before claiming turnkey reproducibility.
