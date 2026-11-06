# Hormonal Modulation of Subjective Sleep Quality and Perivascular Space Volume Across the Menstrual Cycle — Analysis Code

## Repository structure

```
.
├── preprocessing/
│   └── run_fmriprep.sh        # fMRIPrep, anatomical T1w
├── PVSV_calculation/
│   ├── mgz_to_mni.py          # FreeSurfer aparc+aseg -> MNI space
│   ├── run_shiva.sh           # SHiVAi PVS segmentation
│   └── calculate_pvsv.py      # Regional and total PVS volumes from SHiVAi output
├── analysis/
│   └── stats.R                # Statistical analysis
├── plots/
│   ├── fig_style.R            # Shared figure style and helpers
│   ├── fig2_fig4.R            # Fig. 2 (hormone regressions), Fig. 4 (daily dynamics)
│   ├── fig3_moderation.R      # Fig. 3 (gonadotropin moderation, Johnson–Neyman)
│   └── supp_figures.R         # Supplementary Figs. 2–5
├── README.md
└── environment.md             # Software and package versions
```

## Data

[OpenNeuro ds002674](https://openneuro.org/datasets/ds002674/versions/1.0.6)
("28andMe" dense-sampling study). Item-level 1-dPSQI responses, used to
derive the component scores (`psqi_scored`), are not part of the public
release and are available from the dataset's authors on request.

## Imaging pipeline

1. `preprocessing/run_fmriprep.sh`: fMRIPrep 25.1.3, anatomical (T1w)
   preprocessing only, run session-wise via a BIDS filter file.
2. `PVSV_calculation/`: PVS segmentation with
   [SHiVAi](https://github.com/pboutinaud/SHiVAi) 0.4.2 (SHIVA-PVS model)
   and regional volume calculation, in the order `mgz_to_mni.py` →
   `run_shiva.sh` → `calculate_pvsv.py`.
   The SHiVAi container, model weights, configuration file and FreeSurfer
   lookup table referenced in `run_shiva.sh` are available from the
   [SHiVAi repository](https://github.com/pboutinaud/SHiVAi).

## Statistical analysis

`analysis/stats.R` expects a data frame `df` with one row per study day:

- `session`, `condition` (`"natural"` or `"OC"`)
- `total_pvs_volume`, `pvs_volume_centrumsemiovale`, `pvs_volume_basalganglia`
- `psqi` (1-dPSQI total score)
- `` `Estradiol (LCMS; pg/mL)` ``, `` `Progesterone (LCMS; ng/mL)` ``,
  `` `LH (mIU/ml)` ``, `` `FSH (mIU/ml)` ``

and optionally `psqi_scored` (columns `condition`, `session_num`,
`C1_quality` … `C7_daytime`) for the 1-dPSQI components.

The script covers outlier exclusion, condition comparisons (block
permutation, block-bootstrap CIs), hormone regression and interaction
models, the day-to-day association of sleep quality and PVSV (including
regional, lagged and range-restriction analyses) and all sensitivity
analyses (GLS with AR(1) errors, smooth trend of study day, leave-one-out,
rank-based GLS, block lengths 2–5 days). Results are printed as tables;
Q–Q plots and the robustness results are written to the working directory
(`qq_all_models.png`, `autocorrelation_series.csv`,
`robustness_regression.csv`, `robustness_day_to_day.csv`,
`robustness_block_length.csv`).

## Figures

The figure scripts recompute the statistics shown in the panels and write
each figure as vector PDF and 600 dpi TIFF. Supplementary Fig. 1 (example
segmentation) is made from the imaging data.

## Reproducing the analysis

From the repository root, with `df` (and `psqi_scored`) loaded:

```r
source("analysis/stats.R")
source("plots/fig2_fig4.R")
source("plots/fig3_moderation.R")
source("plots/supp_figures.R")
```

## Contact

Daryna Apter — darynaapter@gmail.com
