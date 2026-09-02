# Computational environment

As reported in the manuscript, all analyses were run under:

- **R version 4.5.3** (2026-03-11)

## Required packages

### `analysis/stats.R`

| Package   | Version | Purpose                                             |
|-----------|---------|------------------------------------------------------|
| dplyr     | 1.2.1   | Data wrangling                                        |
| nlme      | 3.1.168 | GLS models with AR(1) correlation structure (`gls()`) |
| lmtest    | 0.9.40  | Durbin–Watson test (`dwtest()`)                        |
| ggplot2   | 4.0.3   | Q–Q plots                                              |
| ggthemes  | 5.2.0   | Plot theme (`theme_foundation()`)                     |
| cowplot   | 1.2.0   | Combining Q–Q plots into a grid (`plot_grid()`)        |

### `plots/delta_reg_plot.R` and `plots/interaction_plot.R`

Additionally to the packages above:

| Package   | Version | Purpose                                             |
|-----------|---------|--------------------------------------------------------|
| gridExtra | 2.3     | Grid-based plot layout (`delta_reg_plot.R` only)        |
| grid      | base R  | Low-level graphics primitives (base package, no install needed) |
| scales    | 1.4.0   | Axis breaks (`pretty_breaks()`), alpha blending          |
| showtext  | 0.9.8   | Unicode-safe font rendering in PDF export                |
| ggtext    | 0.1.2   | Rich-text plot annotations (italics, subscripts)          |

Install with:

```r
install.packages(c("dplyr", "nlme", "lmtest", "ggplot2", "ggthemes", "cowplot",
                    "gridExtra", "scales", "showtext", "ggtext"))
```

