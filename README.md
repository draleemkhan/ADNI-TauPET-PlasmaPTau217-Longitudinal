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

The repository includes **44 MATLAB (`.m`) scripts** and **2 R (`.R`) scripts**. `Source_SHA256.csv` contains **46 script records**, with fingerprints for the historical source scripts and their GitHub copies. The GitHub-copy hashes were updated after directory-path corrections; the original source hashes were preserved. Hash verification establishes file identity, **not** independent execution or reproduction of study results.

## Data access

ADNI source data must be obtained separately by an authorised researcher under the applicable ADNI data-access conditions. Consult [`ADNI_Source_Data_Inventory_14.csv`](01_Documentation/Data_Access/ADNI_Source_Data_Inventory_14.csv) for the documented 14-file source inventory. **Do not upload ADNI source files, participant-level derived datasets, or locally generated participant-level outputs to a public repository.**

**Tau PET source clarification:** For reproduction, use the verified **non-partial-volume-corrected** UC Berkeley 6-mm flortaucipir tau PET source, `UCBERKELEY_TAU_6MM_11Aug2026.csv`. The detailed guide documents a discrepancy between the manuscript's cohort-selection and biomarker-source wording that should be resolved before manuscript publication.

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

Expected **cohort-flow checkpoints** are **1,246 → 1,237 → 1,100 → 469 → 275 → 211 → 210 → 143**; these are validation checkpoints, not counts to hard-code. The final model-ready dataset should contain **143 participants and 379 CDR-SB observations**. Report any deviations rather than adjusting results to match the checkpoints.

## Reproducibility and limitations

- The public repository does **not** contain the restricted ADNI source datasets or participant-level intermediate files needed for execution.
- Some scripts require outputs from preceding scripts. The collection is **not** a single parameterised pipeline and has **not** been independently verified by a complete clean-environment rerun.
- The source-code manifest documents script fingerprints. It does not verify statistical outputs, data permissions, or successful execution on another system.
- The study's reported results should be checked against newly generated outputs; differences should be documented transparently.

For detailed instructions and known caveats, consult the [pipeline reproduction guide](01_Documentation/Analysis_Workflow/ADNI_Pipeline_Reproduction_Guide_Final_46.md).

## Data-use acknowledgement

This project uses data from the **Alzheimer's Disease Neuroimaging Initiative (ADNI)**. Researchers reusing the scripts are responsible for complying with ADNI data-use, acknowledgement, and publication requirements applicable to their own work.
