# Hormonal Modulation of Subjective Sleep Quality and Perivascular Space Volume Across the Menstrual Cycle: A Dense-Sampling Study — Analysis Code

## Repository structure

```
.
├── preprocessing/
│   └── run_fmriprep.sh          # fMRIPrep, anatomical T1w
├── PVSV_calculation/
│   ├── mgz_to_mni.py             # FreeSurfer aparc+aseg -> MNI space
│   ├── run_shiva.sh             # SHiVAi PVS segmentation
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

## Data

[OpenNeuro ds002674](https://openneuro.org/datasets/ds002674/versions/1.0.6)
("28andMe" dense-sampling study). Raw, item-level PSQI responses (used to
derive `psqi_scored`) are not part of the public release; available from
the dataset's authors on request.

## `preprocessing/`

fMRIPrep, anatomical (T1w) preprocessing only, incl. FreeSurfer `aparc+aseg`.
v25.1.3. 

## `PVSV_calculation/`

PVS segmentation ([ShiVAi](https://github.com/pboutinaud/SHiVAi) v0.4.2,
SHIVA-PVS model) and volume calculation from the fMRIPrep/FreeSurfer
output, in pipeline order (`mgz_to_mni.py` → `run_shiva.sh` →
`calculate_pvsv.py`).

## `analysis/stats.R`

Statistical analysis: hormone regressions, PSQI x hormone interaction
models, and day-to-day sleep-PVS coupling (block permutation/bootstrap
throughout). Expects a data frame `df` with one row per study day:

- `session`, `condition` (`"natural"` or `"OC"`)
- `total_pvs_volume`, `pvs_volume_centrumsemiovale`, `pvs_volume_basalganglia`
- `psqi`, optionally `psqi_scored` (`C1_quality` … `C7_daytime`)
- `` `Estradiol (LCMS; pg/mL)` ``, `` `Progesterone (LCMS; ng/mL)` ``,
  `` `LH (mIU/ml)` ``, `` `FSH (mIU/ml)` ``

## `plots/`

Reproduces the manuscript figures from the same data/methods as
`analysis/stats.R`; reuses its results objects if run first in the same R
session, otherwise recomputes locally.

## Reproducing the statistical analysis

```r
# install.packages(c("dplyr", "nlme", "lmtest", "ggplot2", "ggthemes",
#                     "cowplot", "gridExtra", "showtext", "ggtext", "scales"))
source("analysis/stats.R")
source("plots/delta_reg_plot.R")
source("plots/interaction_plot.R")
```

## Contact

Daryna Apter — darynaapter@gmail.com
