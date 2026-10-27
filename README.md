# Overlapping Longitudinal Information From Tau PET and Plasma p-Tau217 for Clinical Progression in Alzheimer's Disease

## Overview

This repository contains MATLAB and R analysis scripts and documentation for a longitudinal study of tau PET, plasma phosphorylated tau 217 (p-tau217), amyloid PET, and clinical progression in amyloid-positive participants from the Alzheimer's Disease Neuroimaging Initiative (ADNI).

The final analytical cohort comprises **143 participants** (33 cognitively normal, 75 with mild cognitive impairment, and 35 with dementia), contributing **379 Clinical Dementia Rating–Sum of Boxes (CDR-SB) observations**. The analyses examine biomarker associations with longitudinal CDR-SB change, incremental model information, and robustness across sensitivity analyses.

**Important:** This is a code-and-documentation repository, **not** a public release of ADNI participant-level data. The scripts have been curated for sharing and their GitHub-copy fingerprints recorded, but an independent, clean-environment, end-to-end rerun has **not** been demonstrated. Do not interpret the repository as a turnkey, one-command pipeline.

## Repository contents

```text
.
├── README.md
├── 01_Documentation/
│   ├── Analysis_Workflow/
│   │   └── ADNI_Pipeline_Reproduction_Guide_Final_46.md
│   └── Data_Access/
│       └── ADNI_Source_Data_Inventory_14.csv
├── 02_Source_Code_Manifest/
│   └── Source_SHA256.csv
├── 03_MATLAB_Codes/
│   ├── 01_Cohort_Construction/
│   ├── 02_Biomarker_Processing/
│   ├── 03_Primary_Models/
│   ├── 04_Incremental_Models/
│   ├── 05_Prespecified_Sensitivity/
│   ├── 06_PostHoc_Analysis/
│   └── 07_Figures_Tables/
└── 04_R_Codes/
```

The repository includes **47 MATLAB (`.m`) scripts** and **2 R (`.R`) scripts**. `Source_SHA256.csv` contains **49 script records**, with fingerprints for the historical source scripts and their GitHub copies. The GitHub-copy hashes were updated after directory-path corrections; the original source hashes were preserved. Hash verification establishes file identity, **not** independent execution or reproduction of study results.

## Data access

ADNI source data must be obtained separately by an authorised researcher under the applicable ADNI data-access conditions. Consult [`ADNI_Source_Data_Inventory_14.csv`](01_Documentation/Data_Access/ADNI_Source_Data_Inventory_14.csv) for the documented 14-file source inventory. **Do not upload ADNI source files, participant-level derived datasets, or locally generated participant-level outputs to a public repository.**

The MIT License applies to the code in this repository only. ADNI data are not included and remain subject to the ADNI Data Use Agreement.

**Tau PET source clarification:** For reproduction, use the UC Berkeley non-partial-volume-corrected (NPVC) 6-mm flortaucipir tau PET dataset, `UCBERKELEY_TAU_6MM_11Aug2026.csv`. This is the source used for cohort construction and baseline meta-temporal tau PET quantification. The partial-volume-corrected dataset (`UCBERKELEY_TAUPVC_6MM`) was not used.

## Software and local paths

The analysis was developed using **MATLAB R2025a** (including mixed-effects modelling functionality such as `fitlme`) and **R** (including the **geepack** package for GEE analyses). Additional dependencies should be checked in the individual scripts.

The curated GitHub copies use the following Windows reproduction conventions:

| Location | Purpose |
| --- | --- |
| `D:\source files` | Authorised original ADNI source CSVs |
| `D:\codes` | Optional local working location for downloaded scripts |
| `D:\output` | Intermediate inputs and generated outputs, **including TIFF, PDF, and PNG figures directly in this folder** |

The repository itself organises scripts in nested folders under `03_MATLAB_Codes/`. **Do not assume that flattening scripts into `D:\codes` preserves relative-path dependencies.** In particular, `Step11B_Broad_ADNI_Registry_Activity_Audit.m` has a repository-relative lookup for the source-data inventory; its compatibility with a flattened `D:\codes` layout has not been confirmed. Review each script's actual input and output dependencies before execution.

## Reproduction workflow

Start with the **[detailed pipeline reproduction guide](01_Documentation/Analysis_Workflow/ADNI_Pipeline_Reproduction_Guide_Final_46.md)**. It describes the expected input files, cohort construction, biomarker standardisation, modelling specifications, dependency order, robustness analyses, and verification checkpoints.

At a high level, the workflow covers:

1. Construct the index flortaucipir tau PET cohort, match amyloid PET and plasma p-tau217, and apply clinical and longitudinal eligibility criteria.
2. Prepare the final participant-level biomarkers and repeated CDR-SB observations, using time relative to index tau PET.
3. Fit the primary longitudinal mixed-effects model with biomarker-by-time interactions and uncorrelated participant-specific random intercepts and slopes.
4. Compare nested biomarker models using maximum likelihood and evaluate incremental model information.
5. Run the documented sensitivity, post hoc, attrition, and figure/table scripts in their required dependency order.

**Amyloid PET tracer eligibility:** Initial screening permitted [¹⁸F]florbetapir, [¹⁸F]florbetaben, and [¹⁸F]NAV4694. After applying the prespecified imaging, biomarker, clinical, and longitudinal eligibility criteria, the final analytical cohort of 143 participants included 74 florbetapir and 69 florbetaben examinations. No NAV4694 examinations were represented in the final cohort. NAV4694 was not retrospectively excluded from the original eligibility criteria.

Expected **cohort-flow checkpoints** are **1,246 → 1,237 → 1,100 → 469 → 275 → 211 → 210 → 143**; these are validation checkpoints, not counts to hard-code. The final model-ready dataset should contain **143 participants and 379 CDR-SB observations**. Report any deviations rather than adjusting results to match the checkpoints.

**Locked analysis workspace:** Run `Step08C_Construct_Locked_N143_Workspace.m` immediately after Step08. It rebuilds `ADNI_REVISED_COHORT_LOCKED_N143.mat`, which downstream scripts load; the reconstruction was verified identical to the original locked workspace.

## Reproducibility and limitations

- The public repository does **not** contain the restricted ADNI source datasets or participant-level intermediate files needed for execution.
- Some scripts require outputs from preceding scripts. The collection is **not** a single parameterised pipeline and has **not** been independently verified by a complete clean-environment rerun.
- The source-code manifest documents script fingerprints. It does not verify statistical outputs, data permissions, or successful execution on another system.
- The study's reported results should be checked against newly generated outputs; differences should be documented transparently.

For detailed instructions and known caveats, consult the [pipeline reproduction guide](01_Documentation/Analysis_Workflow/ADNI_Pipeline_Reproduction_Guide_Final_46.md).

## Data-use acknowledgement

This project uses data from the **Alzheimer's Disease Neuroimaging Initiative (ADNI)**. Researchers reusing the scripts are responsible for complying with ADNI data-use, acknowledgement, and publication requirements applicable to their own work.
