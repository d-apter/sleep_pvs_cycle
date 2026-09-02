# Hormonal Modulation of Subjective Sleep Quality and Perivascular Space Volume Across the Menstrual Cycle: A Dense-Sampling Study — Analysis Code

## Overview

Analysis of sleep quality (PSQI) and perivascular space volume (PVSV)
across the physiological menstrual cycle and oral contraceptive (OCP) use,
using data from [OpenNeuro ds002674](https://openneuro.org/datasets/ds002674/versions/1.0.6)
("28andMe" dense-sampling study). Covers the full pipeline from T1w images
to statistical results: anatomical preprocessing, PVS segmentation, PVS
volume calculation, and the statistical analysis (`analysis/stats.R`).

Because observations are repeated daily measurements within one
individual, they are not independent across days. The statistical pipeline
therefore relies on **moving-block permutation and bootstrap procedures**
(block length `L = n^(1/3)`) rather than standard parametric tests
wherever the independence assumption would be violated:

1. **Hormone → sleep quality / PVSV** (Step 1): linear models (OLS or GLS
   AR(1)) per hormone (estradiol, progesterone, LH, FSH), physiological
   cycle only.
2. **Hormone × sleep quality interaction models** (Step 2) predicting
   PVSV. Steps 1–2 share diagnostics (DW, Shapiro–Wilk, Q–Q, Cook's
   distance) and a leave-one-out check for FDR-significant models.
3. **Sleep–PVS coupling** (Step 3, day-to-day change level): Spearman
   correlations between ΔPSQI and ΔPVSV per condition, region-specific
   (CSO/BG) decomposition, partial correlations, a difference-in-coupling
   test between conditions, and lagged (t+1 to t+3) checks.

P-values from repeated tests of the same outcome (Steps 1–2) are FDR-corrected
(Benjamini–Hochberg) within each model family.

## Repository structure

```
.
├── preprocessing/
│   └── run_fmriprep.sh          # fMRIPrep, anatomical (T1w) only
├── PVSV_calculation/
│   ├── mgz_to_mni.py             # FreeSurfer aparc+aseg -> MNI space (SHiVAi custom LUT input)
│   ├── run_shiva.sh             # SHiVAi PVS segmentation (Apptainer/Singularity)
│   └── calculate_pvsv.py        # Regional + total PVS volumes from SHiVAi output
├── analysis/
│   └── stats.R                  # Statistical analysis
├── plots/
│   ├── delta_reg_plot.R         # Figure 1 (trajectories + delta coupling) & Figure 2 (hormone regressions)
│   └── interaction_plot.R       # Figure: PSQI x hormone interaction + conditional slopes
├── README.md
├── .gitignore
└── environment.md               # R version and package versions used
```

## Imaging: preprocessing & PVS segmentation

Anatomical preprocessing (incl. FreeSurfer `aparc+aseg` parcellation) was
performed with [fMRIPrep](https://fmriprep.org/) v25.1.3. PVS detection
used a containerised [ShiVAi](https://github.com/pboutinaud/SHiVAi) (v0.4.2,
`shiva_0.4.2.sif`), based on the SHIVA-PVS U-Net model (weights openly
available in the ShiVAi repo).

Pipeline order:

- **`preprocessing/run_fmriprep.sh`**: fMRIPrep, anatomical (T1w) only
  (`--anat-only`).
- **`PVSV_calculation/mgz_to_mni.py`**: warps each session's FreeSurfer
  `aparc+aseg` into MNI space, used as SHiVAi's custom brain segmentation
  (`--custom_LUT`).
- **`PVSV_calculation/run_shiva.sh`**: the SHiVAi segmentation call.
- **`PVSV_calculation/calculate_pvsv.py`**: aggregates SHiVAi's per-region
  output (`pvs_stats.csv`) into `pvs_volume_centrumsemiovale`,
  `pvs_volume_basalganglia`, `pvs_volume_hippocampus`, and
  `total_pvs_volume` (sum of the three).

`run_fmriprep.sh` paths are placeholders (`/path/to/...`) to edit before
use. `mgz_to_mni.py`/`calculate_pvsv.py` assume a small set of
project-specific helper variables/functions (subject list, dataset table,
path helpers, cluster job submission).

The underlying imaging data and SHiVAi output (`pvs_stats.csv`) are not
redistributed here. The resulting `total_pvs_volume`,
`pvs_volume_centrumsemiovale`, and `pvs_volume_basalganglia` values used
as input to `analysis/stats.R` are available from the corresponding
author of the manuscript on reasonable request.

## Statistical analysis data

`analysis/stats.R` expects a data frame `df` with one row per study day:

- `session`, `condition` (`"natural"` or `"OC"`)
- `total_pvs_volume`, `pvs_volume_centrumsemiovale`, `pvs_volume_basalganglia`
- `psqi` (PSQI global score; optionally `psqi_scored` with sub-component
  columns `C1_quality` … `C7_daytime` — not publicly available, see note
  below)
- `` `Estradiol (LCMS; pg/mL)` ``, `` `Progesterone (LCMS; ng/mL)` ``,
  `` `LH (mIU/ml)` ``, `` `FSH (mIU/ml)` ``

Raw data are publicly available from
[OpenNeuro ds002674](https://openneuro.org/datasets/ds002674/versions/1.0.6)
and are not redistributed here. See `environment.md` and the manuscript's
Methods section for details, e.g. PSQI sub-component scoring.

**Note on PSQI item-level data**: `psqi` (global/total score) is in the
public OpenNeuro release. The raw, item-level responses used to derive
`psqi_scored` (Family 3 sub-component comparison) are **not** public;
these were obtained directly from the dataset's authors upon request.

## Figures

`plots/delta_reg_plot.R` and `plots/interaction_plot.R` reproduce the main
manuscript figures using the same methods as `analysis/stats.R`
(permutation/bootstrap and model-fitting logic duplicated locally in each
script). Where possible, both reuse objects already computed by
`analysis/stats.R` (e.g. `res_delta_nat`, `ci_nat`) if that script ran
first in the same R session, so figures match the reported numbers
exactly; otherwise they recompute locally.

- **`delta_reg_plot.R`**: hormone/PSQI/PVSV trajectories, day-to-day
  (delta) sleep-PVS coupling scatterplots, hormone regression scatterplots
  (PVSV ~ Estradiol, PSQI ~ LH).
- **`interaction_plot.R`**: PSQI x hormone interaction plots (PVSV ~ PSQI x
  LH/FSH) and conditional-slope trajectories.

## Reproducing the statistical analysis

```r
# install.packages(c("dplyr", "nlme", "lmtest", "ggplot2", "ggthemes",
#                     "cowplot", "gridExtra", "showtext", "ggtext", "scales"))
source("analysis/stats.R")
source("plots/delta_reg_plot.R")
source("plots/interaction_plot.R")
```

`set.seed(42)` is set globally and reset at the start of every resampling
function, so results are reproducible regardless of run order.

## Contact

Daryna Apter — darynaapter@gmail.com
