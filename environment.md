# Computational environment

All analyses were run under **R version 4.5.3** (2026-03-11).

## R packages

| Package | Version | Used in | Purpose |
|---------|---------|---------|---------|
| dplyr   | 1.2.1   | all     | Data wrangling |
| nlme    | 3.1.168 | `analysis/stats.R` | GLS models with AR(1) errors (`gls()`) |
| lmtest  | 0.9.40  | `analysis/stats.R` | Durbin–Watson test (`dwtest()`) |
| ggplot2 | 4.0.3   | all     | Q–Q plots and figures |
| cowplot | 1.2.0   | all     | Panel layout |
| scales  | 1.4.0   | `plots/` | Axis breaks and labels |

`grid`, `grDevices` and `stats` are part of base R. Figures use the Arial
font (embedded via the Quartz device on macOS or Cairo elsewhere).

```r
install.packages(c("dplyr", "nlme", "mgcv", "lmtest", "ggplot2", "cowplot", "scales"))
```

## Imaging and Python

| Software | Version | Purpose |
|----------|---------|---------|
| fMRIPrep (Singularity image) | 25.1.3 | Anatomical preprocessing, FreeSurfer `aparc+aseg` |
| SHiVAi (Singularity image) | 0.4.2 | PVS segmentation (SHIVA-PVS model) |
| Python with pandas | – | Aggregation of regional PVS volumes |
