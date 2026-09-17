# Grip strength and falls in four national cohorts: analysis code

This repository contains the analysis code for the manuscript

> **Grip strength and falls in 114,129 middle-aged and older women across four national cohorts:
> consistent observational associations but no replication in an independent Mendelian
> randomisation analysis**

## What the code does

The pipeline has four parts.

| Folder | Contents |
| --- | --- |
| `01_build` | Builds the CHARLS analysis dataset from the harmonised CHARLS release and the raw 2011–2018 wave files; also contains the shared `db_paths.R` used to locate the source databases |
| `02_cohort` | The four-cohort cohort analysis: sample assembly for CHARLS, ELSA, HRS and SHARE, design-based survey models, DerSimonian–Laird pooling, design effects, Table 1, Table 2 and subgroups, dose–response splines, competing-risk analyses and all cohort figures |
| `03_mendelian_randomisation` | Two-sample MR: the overlapping versus independent outcome contrast, sensitivity estimators, instrument restriction, positive and negative controls, the reverse-direction analysis and Figure 4 |
| `04_supplementary_analyses` | Analyses added at revision: the CHARLS age-at-natural-menopause effect-modification analysis and the first-stage sampling unit sensitivity analysis |
| `05_quality_checks` | `check_manuscript_numbers.py` re-derives every key number in the manuscript from the result files and exits non-zero if any value disagrees; `build_report.py` writes the audit table into the report |

## Data

None of the source data are redistributed here.

* **Cohort data.** CHARLS, ELSA, HRS and SHARE are obtained from their own data providers under
  data use agreements; harmonised versions come from the Gateway to Global Aging Data. The scripts
  expect them under the paths defined in `01_build/db_paths.R`.
* **OpenGWAS access.** The reverse-direction scripts read an OpenGWAS JSON web token from
  `/tmp/opengwas_jwt`; obtain your own token from OpenGWAS and write it to that path (or edit the
  path in those five scripts). No credentials are stored in this repository.
* **Genetic data.** All summary statistics are publicly available: IEU OpenGWAS accessions
  `ukb-b-7478`, `ukb-b-10215`, `ebi-a-GCST90007526`, `ukb-b-2535`, `ebi-a-GCST90012857`,
  `ieu-b-40`, `ieu-a-300`, and the FinnGen release 12 endpoint `FALLS`.

## Environment

R 4.6.1 with the following packages (versions used for the submitted manuscript):

```
survey        4.5
data.table    1.18.6.1
TwoSampleMR   0.7.9
ieugwasr      1.1.0
MRPRESSO      1.0
haven         2.5.5
ggplot2       4.0.3
ragg          1.5.2
dplyr         1.2.1
gridExtra     2.3.1
```

Python 3.12 is needed only for the quality checks (`05_quality_checks`), which use the standard
library alone.

## How to run

The scripts were written as a research pipeline rather than as a package, and they use **absolute
paths** pointing at the author's local environment. To reproduce the analysis elsewhere:

1. Edit the path constants at the top of each script (they are declared as `B`, `OUT`, `RES`, `W` or
   similar) to point at your own copies of the data and at a writable output directory.
2. Edit `01_build/db_paths.R` so that it resolves your copies of the four cohort databases.
3. Run the folders in numerical order: `01_build` → `02_cohort` (in file-name order; `00_deff_helper.R`
   must be sourced before the cohort scripts, which they do themselves) → `03_mendelian_randomisation`
   → `04_supplementary_analyses`.
4. After the pipeline has run, execute `05_quality_checks/check_manuscript_numbers.py` against the
   manuscript folder; it re-computes and compares every reported number.

Two differences from a conventional pipeline are worth flagging. First, the figures are drawn with
`ragg`, and the figure scripts accept an environment variable `FIGLANG` (`FIGLANG=en`) that switches
all labels to English. Second, `02_cohort/10_charls_design_effect.R` both fits the CHARLS primary
model and exports the analysis sample (`t1/样本_CHARLS.rds`) that the other cohort scripts and the
subgroup, dose–response and covariate-specification scripts read.

## Note on the analysis history

`04_supplementary_analyses` and the supplementary analysis plan document the changes made after the
original analysis plan: the ELSA clustering unit was revised during an internal methodological
audit, the first-stage sampling unit for CHARLS was re-examined, the age-at-menopause analysis was
added, and a within-person grip-strength slope covariate was removed from the CHARLS primary model so
that all four cohorts share one analytical template. The manuscript reports each of these changes
explicitly.
