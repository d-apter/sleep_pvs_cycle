# Title:  Statistical analysis - sleep quality, hormones, and PVSV
# Author: Daryna Apter
# Date:   2026-08-31

# Analysis pipeline: outlier check, descriptives, condition comparison,
# hormone regression (Step 1), moderation models (Step 2), diagnostics,
# sleep-PVS delta coupling (Step 3).

library(dplyr)
library(stats)
filter <- dplyr::filter
select <- dplyr::select
arrange <- dplyr::arrange
mutate <- dplyr::mutate
group_by <- dplyr::group_by
summarise <- dplyr::summarise
library(nlme)
library(lmtest)
library(ggplot2)
library(ggthemes)
library(cowplot)

# setwd("/path/to/project")  # edit as needed

# df must contain the columns described in README.md, with condition
# values "natural"/"OC"
# df <- read.csv("data.csv")

set.seed(42)

# ---- 0. DATA PREPARATION ----

# Data-level condition values are "natural"/"OC"; "OCP" is used in labels only.
nat <- df %>%
  dplyr::filter(condition == "natural") %>%
  arrange(session) %>%
  mutate(time = seq_len(n()))

oc <- df %>%
  dplyr::filter(condition == "OC") %>%
  arrange(session) %>%
  mutate(time = seq_len(n()))

# Log-transformation of hormones (all steroids + gonadotropins)
nat <- nat %>%
  mutate(
    log_E2  = log(`Estradiol (LCMS; pg/mL)`   + 1),
    log_P4  = log(`Progesterone (LCMS; ng/mL)` + 1),
    log_LH  = log(`LH (mIU/ml)`  + 1),
    log_FSH = log(`FSH (mIU/ml)` + 1)
  )
oc <- oc %>%
  mutate(
    log_E2  = log(`Estradiol (LCMS; pg/mL)`   + 1),
    log_P4  = log(`Progesterone (LCMS; ng/mL)` + 1),
    log_LH  = log(`LH (mIU/ml)`  + 1),
    log_FSH = log(`FSH (mIU/ml)` + 1)
  )

# Block sizes for permutation/bootstrap (rule of thumb b ~ n^(1/3))
DEFAULT_BLOCK_SIZE   <- 3
DIFF_TEST_BLOCK_SIZE <- 4

# ---- 1. OUTLIER CHECK (|z| > 3, pre-registered criterion) ----

# Outlier check across raw and delta values of PVSV/PSQI, both conditions
check_outlier <- function(x) sum(abs(scale(x)[,1]) > 3, na.rm = TRUE)

n_out_raw <- c(
  nat_pvs  = check_outlier(nat$total_pvs_volume),
  nat_psqi = check_outlier(nat$psqi),
  oc_pvs   = check_outlier(oc$total_pvs_volume),
  oc_psqi  = check_outlier(oc$psqi)
)
n_out_delta <- c(
  nat_dpvs  = check_outlier(c(NA, diff(nat$total_pvs_volume))),
  nat_dpsqi = check_outlier(c(NA, diff(nat$psqi))),
  oc_dpvs   = check_outlier(c(NA, diff(oc$total_pvs_volume))),
  oc_dpsqi  = check_outlier(c(NA, diff(oc$psqi)))
)

n_out_raw
n_out_delta

# Pre-registered |z|>3 outlier removed from all analyses; delta analyses
# (Step 3) also drop the following day.
z_psqi_oc       <- scale(oc$psqi)[, 1]
outlier_session <- oc$session[which(abs(z_psqi_oc) > 3)]

oc_clean <- oc %>% filter(!(session %in% outlier_session))

# Significance stars, used throughout
sig_star <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01,  "**",
                ifelse(p < 0.05,  "*",
                       ifelse(p < 0.10,  "\u2020", "ns"))))
}

# ---- 2. DESCRIPTIVE STATISTICS (Table 1) ----
desc_vars <- list(
  list(v = "total_pvs_volume",           lab = "PVS volume, total (mm3)"),
  list(v = "pvs_volume_centrumsemiovale", lab = "PVS volume, CSO (mm3)"),
  list(v = "pvs_volume_basalganglia",     lab = "PVS volume, BG (mm3)"),
  list(v = "psqi",                        lab = "PSQI global score"),
  list(v = "Estradiol (LCMS; pg/mL)",    lab = "Estradiol (pg/mL)"),
  list(v = "Progesterone (LCMS; ng/mL)", lab = "Progesterone (ng/mL)"),
  list(v = "LH (mIU/ml)",                lab = "LH (mIU/mL)"),
  list(v = "FSH (mIU/ml)",               lab = "FSH (mIU/mL)")
)

descriptive_stats <- do.call(rbind, lapply(desc_vars, function(v) {
  x_nat <- nat[[v$v]]
  x_oc  <- oc_clean[[v$v]]
  data.frame(
    variable   = v$lab,
    nat_mean   = mean(x_nat, na.rm = TRUE), nat_sd  = sd(x_nat, na.rm = TRUE),
    nat_min    = min(x_nat,  na.rm = TRUE), nat_max = max(x_nat, na.rm = TRUE),
    oc_mean    = mean(x_oc,  na.rm = TRUE), oc_sd   = sd(x_oc,  na.rm = TRUE),
    oc_min     = min(x_oc,   na.rm = TRUE), oc_max  = max(x_oc,  na.rm = TRUE)
  )
}))

descriptive_stats

# ---- 3. CONDITION COMPARISON (permutation, labels across all 60 sessions) ----

# Combined dataset (n = 59: 30 physiological cycle + 29 OCP, outlier removed)
df_combined <- bind_rows(
  nat       %>% mutate(condition = "natural"),
  oc_clean  %>% mutate(condition = "OC")
)

# Block-bootstrap CI for the independent two-group mean difference
block_bootstrap_ci_meandiff <- function(x1, x2,
                                        block_size1 = DEFAULT_BLOCK_SIZE,
                                        block_size2 = DEFAULT_BLOCK_SIZE,
                                        n_boot = 10000, conf = 0.95) {
  set.seed(42)  # reproducibility lock (see run_condition_perm)
  n1 <- length(x1); n2 <- length(x2)
  block_id1 <- ceiling(seq_len(n1) / block_size1); n_blocks1 <- max(block_id1)
  block_id2 <- ceiling(seq_len(n2) / block_size2); n_blocks2 <- max(block_id2)
  
  diffs_boot <- replicate(n_boot, {
    boot_blocks1 <- sample(seq_len(n_blocks1), size = n_blocks1, replace = TRUE)
    boot_idx1    <- unlist(lapply(boot_blocks1, function(b) which(block_id1 == b)))[seq_len(n1)]
    boot_blocks2 <- sample(seq_len(n_blocks2), size = n_blocks2, replace = TRUE)
    boot_idx2    <- unlist(lapply(boot_blocks2, function(b) which(block_id2 == b)))[seq_len(n2)]
    mean(x1[boot_idx1], na.rm = TRUE) - mean(x2[boot_idx2], na.rm = TRUE)
  })
  alpha <- 1 - conf
  c(lo = as.numeric(quantile(diffs_boot, alpha / 2)),
    hi = as.numeric(quantile(diffs_boot, 1 - alpha / 2)))
}

# Moving-block permutation of the mean difference
run_condition_perm <- function(var, label, block_size = DIFF_TEST_BLOCK_SIZE,
                               n_perm = 10000) {
  set.seed(42)  # reset every call for order-independent reproducibility
  x <- df_combined[[var]]
  g <- df_combined$condition
  n_nat <- sum(g == "natural")
  n_total <- length(x)
  
  obs_diff <- mean(x[g == "natural"], na.rm = TRUE) - mean(x[g == "OC"], na.rm = TRUE)
  
  block_id <- ceiling(seq_len(n_total) / block_size)
  n_blocks <- max(block_id)
  
  perm_diffs <- replicate(n_perm, {
    block_order   <- sample(seq_len(n_blocks))
    reordered_idx <- unlist(lapply(block_order, function(b) which(block_id == b)))
    x_reordered   <- x[reordered_idx]
    mean(x_reordered[1:n_nat], na.rm = TRUE) -
      mean(x_reordered[(n_nat + 1):n_total], na.rm = TRUE)
  })
  p_perm <- (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1)
  
  ci <- block_bootstrap_ci_meandiff(x[g == "natural"], x[g == "OC"],
                                    block_size1 = DEFAULT_BLOCK_SIZE,
                                    block_size2 = DEFAULT_BLOCK_SIZE)
  
  list(var = var, label = label, diff = obs_diff,
       ci_lo = ci["lo"], ci_hi = ci["hi"], p_perm = p_perm)
}

cond_vars_hormone <- list(
  list(v = "Estradiol (LCMS; pg/mL)",    lab = "Estradiol"),
  list(v = "Progesterone (LCMS; ng/mL)", lab = "Progesterone"),
  list(v = "LH (mIU/ml)",                lab = "LH"),
  list(v = "FSH (mIU/ml)",               lab = "FSH")
)

cond_vars_outcome <- list(
  list(v = "total_pvs_volume",           lab = "PVS volume (total)"),
  list(v = "pvs_volume_centrumsemiovale", lab = "PVS volume (CSO)"),
  list(v = "pvs_volume_basalganglia",     lab = "PVS volume (BG)"),
  list(v = "psqi",                        lab = "PSQI global score")
)

# No FDR correction: each variable tests an independent question
cond_results_hormone <- lapply(cond_vars_hormone, function(v) run_condition_perm(v$v, v$lab))
cond_results_outcome <- lapply(cond_vars_outcome, function(v) run_condition_perm(v$v, v$lab))
cond_results_hormone
cond_results_outcome

# Family 3: PSQI sub-components. Requires `psqi_scored`. C6 excluded
# (zero variance); PSQI_total not retested (already in Family 2).
run_condition_perm_generic <- function(data, var, label,
                                       condition_col = "condition",
                                       ref_level = "natural",
                                       block_size = DIFF_TEST_BLOCK_SIZE,
                                       n_perm = 10000) {
  set.seed(42)  # reproducibility lock (see run_condition_perm)
  x <- data[[var]]
  g <- data[[condition_col]]
  n_ref   <- sum(g == ref_level)
  n_total <- length(x)
  
  obs_diff <- mean(x[g == ref_level], na.rm = TRUE) -
    mean(x[g != ref_level], na.rm = TRUE)
  
  block_id <- ceiling(seq_len(n_total) / block_size)
  n_blocks <- max(block_id)
  
  perm_diffs <- replicate(n_perm, {
    block_order   <- sample(seq_len(n_blocks))
    reordered_idx <- unlist(lapply(block_order, function(b) which(block_id == b)))
    x_reordered   <- x[reordered_idx]
    mean(x_reordered[1:n_ref], na.rm = TRUE) -
      mean(x_reordered[(n_ref + 1):n_total], na.rm = TRUE)
  })
  p_perm <- (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1)
  
  list(var = var, label = label, diff = obs_diff, p_perm = p_perm)
}

if (exists("psqi_scored")) {
  psqi_scored_clean <- psqi_scored %>%
    filter(!(condition == "OC" & session_num %in% outlier_session))
  
  psqi_scored_ordered <- psqi_scored_clean %>%
    arrange(factor(condition, levels = c("natural", "OC")), session_num)
  
  subcomp_vars <- list(
    list(v = "C1_quality",      lab = "C1 Subjective quality"),
    list(v = "C2_latency",      lab = "C2 Sleep latency"),
    list(v = "C3_duration",     lab = "C3 Sleep duration"),
    list(v = "C4_efficiency",   lab = "C4 Sleep efficiency"),
    list(v = "C5_disturbances", lab = "C5 Sleep disturbances"),
    list(v = "C7_daytime",      lab = "C7 Daytime dysfunction")
  )
  
  cond_results_subcomp <- lapply(subcomp_vars, function(v)
    run_condition_perm_generic(psqi_scored_ordered, v$v, v$lab))
  cond_results_subcomp
}

# ---- HELPER FUNCTIONS ----

# GLS AR(1) or OLS, chosen via the Durbin-Watson test.
fit_ar <- function(formula_str, data) {
  m0 <- lm(as.formula(formula_str), data = data)
  dw <- dwtest(m0, alternative = "two.sided")
  if (dw$p.value < 0.05) {
    m <- gls(as.formula(formula_str),
             data        = data,
             correlation = corAR1(form = ~time),
             method      = "ML")
    list(model  = m,
         type   = "GLS AR(1)",
         dw_p   = dw$p.value,
         dw_stat= dw$statistic,
         resid  = residuals(m, type = "normalized"))  # AR(1)-adjusted
  } else {
    list(model  = m0,
         type   = "OLS",
         dw_p   = dw$p.value,
         dw_stat= dw$statistic,
         resid  = residuals(m0))
  }
}

# Extract coefficients + CI (works for both GLS and OLS)
get_coef <- function(fit, term) {
  m <- fit$model
  if (fit$type == "GLS AR(1)") {
    s   <- summary(m)$tTable
    b   <- s[term, "Value"]
    se  <- s[term, "Std.Error"]
    p   <- s[term, "p-value"]
    df_ <- s[term, "DF"]
    ci  <- b + c(-1,1) * qt(0.975, df_) * se
  } else {
    s   <- summary(m)$coefficients
    b   <- s[term, "Estimate"]
    se  <- s[term, "Std. Error"]
    p   <- s[term, "Pr(>|t|)"]
    ci  <- confint(m)[term,]
  }
  list(beta = b, se = se, p = p, ci_lo = ci[1], ci_hi = ci[2])
}

# Diagnostics for one model
diag_model <- function(fit, label) {
  res <- fit$resid
  sw  <- shapiro.test(res)
  list(label    = label,
       type     = fit$type,
       dw_stat  = fit$dw_stat,
       dw_p     = fit$dw_p,
       sw_W     = sw$statistic,
       sw_p     = sw$p.value)
}

hormones  <- c("log_E2","log_P4","log_LH","log_FSH")
horm_labs <- c("Estradiol","Progesterone","LH","FSH")

# ---- STEP 1: HORMONE -> PSQI & HORMONE -> PVS ----

# 1a: Hormone -> PSQI
fits_psqi <- list()
diag_psqi <- list()

for (i in seq_along(hormones)) {
  fml  <- paste0("psqi ~ ", hormones[i])
  fit  <- fit_ar(fml, nat)
  fits_psqi[[hormones[i]]] <- fit
  diag_psqi[[i]] <- diag_model(fit, paste0("PSQI ~ ", horm_labs[i]))
}

# 1b: Hormone -> PVS
fits_pvs <- list()
diag_pvs <- list()

for (i in seq_along(hormones)) {
  fml  <- paste0("total_pvs_volume ~ ", hormones[i])
  fit  <- fit_ar(fml, nat)
  fits_pvs[[hormones[i]]] <- fit
  diag_pvs[[i]] <- diag_model(fit, paste0("PVS ~ ", horm_labs[i]))
}

# FDR correction, Step 1
p_psqi <- sapply(hormones, function(h) {
  get_coef(fits_psqi[[h]], h)$p
})
p_pvs <- sapply(hormones, function(h) {
  get_coef(fits_pvs[[h]], h)$p
})

p_psqi_fdr <- p.adjust(p_psqi, method = "BH")
p_pvs_fdr  <- p.adjust(p_pvs,  method = "BH")
p_psqi
p_pvs
p_psqi_fdr
p_pvs_fdr

# ---- STEP 2: MODERATION, PVS ~ PSQI x HORMONE ----
fits_mod  <- list()
diag_mod  <- list()
p_inter   <- numeric(length(hormones))
names(p_inter) <- hormones

# Mean-centre before forming the interaction term (reduces collinearity)
nat_c <- nat %>% mutate(psqi_c = psqi - mean(psqi))
for (h in hormones) {
  nat_c[[paste0(h, "_c")]] <- nat_c[[h]] - mean(nat_c[[h]])
}

for (i in seq_along(hormones)) {
  h_c <- paste0(hormones[i], "_c")
  fml <- paste0("total_pvs_volume ~ psqi_c * ", h_c)
  fit <- fit_ar(fml, nat_c)
  fits_mod[[hormones[i]]] <- fit
  diag_mod[[i]] <- diag_model(fit,
                              paste0("PVS ~ PSQIx", horm_labs[i]))
  
  inter_term <- paste0("psqi_c:", h_c)
  cf <- get_coef(fit, inter_term)
  p_inter[i] <- cf$p
}

# FDR across all interaction terms
p_inter_fdr <- p.adjust(p_inter, method = "BH")
p_inter
p_inter_fdr

# Cook's distance (threshold 4/n, OLS models only)
cooks_diag <- function(fit, label, coef_term, cooks_thresh_n) {
  if (fit$type != "OLS") {
    return(list(label = label, applicable = FALSE))
  }
  cd <- cooks.distance(fit$model)
  thresh <- 4 / cooks_thresh_n
  list(label = label, applicable = TRUE,
       max_cd = max(cd), max_idx = which.max(cd),
       thresh = thresh, n_influential = sum(cd > thresh), n = length(cd))
}

cooks_distance_diag <- list(
  psqi_hormone = lapply(seq_along(hormones), function(i)
    cooks_diag(fits_psqi[[hormones[i]]], paste0("PSQI ~ ", horm_labs[i]),
               hormones[i], nrow(nat))),
  pvs_hormone = lapply(seq_along(hormones), function(i)
    cooks_diag(fits_pvs[[hormones[i]]], paste0("PVS ~ ", horm_labs[i]),
               hormones[i], nrow(nat))),
  interaction = lapply(seq_along(hormones), function(i) {
    h_c <- paste0(hormones[i], "_c")
    cooks_diag(fits_mod[[hormones[i]]], paste0("PVS ~ PSQIx", horm_labs[i]),
               paste0("psqi_c:", h_c), nrow(nat_c))
  })
)

cooks_distance_diag

# Leave-one-out check for FDR-significant hormone models
sig_hormones_psqi <- hormones[p_psqi_fdr < 0.05]
if (length(sig_hormones_psqi) == 0) sig_hormones_psqi <- hormones[p_psqi < 0.05]
sig_hormones_pvs  <- hormones[p_pvs_fdr < 0.05]

run_loo_check <- function(fit_full, formula_str, data, term, label) {
  cd     <- cooks.distance(fit_full$model)
  thresh <- 4 / nrow(data)
  max_idx <- which.max(cd)
  max_cd  <- cd[max_idx]
  influential_session <- data$session[max_idx]
  
  if (max_cd <= thresh) {
    return(list(label = label, loo_run = FALSE,
                max_cd = max_cd, thresh = thresh))
  }
  
  day_excl <- influential_session
  data_loo <- data %>% filter(session != day_excl) %>% mutate(time = seq_len(n()))
  fit_loo <- fit_ar(formula_str, data_loo)
  cf_full <- get_coef(fit_full, term)
  cf_loo  <- get_coef(fit_loo, term)
  
  list(label = label, loo_run = TRUE, day_excluded = day_excl,
       max_cd = max_cd, thresh = thresh,
       beta_full = cf_full$beta, p_full = cf_full$p,
       beta_loo = cf_loo$beta, se_loo = cf_loo$se, p_loo = cf_loo$p,
       ci_lo_loo = cf_loo$ci_lo, ci_hi_loo = cf_loo$ci_hi,
       n_full = nrow(data), n_loo = nrow(data_loo))
}

loo_results <- list()
for (h in intersect(hormones, sig_hormones_psqi)) {
  i <- which(hormones == h)
  loo_results[[paste0("psqi_", h)]] <- run_loo_check(
    fits_psqi[[h]], paste0("psqi ~ ", h), nat, h,
    label = paste0("PSQI ~ ", horm_labs[i]))
}
for (h in intersect(hormones, sig_hormones_pvs)) {
  i <- which(hormones == h)
  loo_results[[paste0("pvs_", h)]] <- run_loo_check(
    fits_pvs[[h]], paste0("total_pvs_volume ~ ", h), nat, h,
    label = paste0("PVS ~ ", horm_labs[i]))
}

loo_results

# No LOO check for Step 2 interaction terms (would contradict moderation logic)

# ---- Q-Q PLOTS FOR ALL REGRESSION MODELS ----
make_qq <- function(res, label, type) {
  df_qq <- data.frame(
    theoretical = qqnorm(res, plot.it = FALSE)$x,
    sample      = qqnorm(res, plot.it = FALSE)$y
  )
  sw  <- shapiro.test(res)
  
  ggplot(df_qq, aes(x = theoretical, y = sample)) +
    geom_abline(slope = 1, intercept = 0,
                color = "grey60", linewidth = 0.5, linetype = "dashed") +
    geom_point(size = 1.8, shape = 21,
               fill = "#3A8C8C", color = "white",
               stroke = 0.3, alpha = 0.85) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("W=%.3f, p=%.3f %s\n%s",
                             sw$statistic, sw$p.value,
                             sig_star(sw$p.value), type),
             hjust = -0.05, vjust = 1.3,
             size = 2.5, color = "grey20", lineheight = 1.2) +
    labs(title = label,
         x = "Theoretical quantiles",
         y = "Sample quantiles") +
    theme_foundation(base_size = 8, base_family = "sans") +
    theme(
      plot.title       = element_text(face = "bold", size = 7.5,
                                      hjust = 0, color = "grey10"),
      axis.title       = element_text(size = 7.5, color = "grey20"),
      axis.text        = element_text(size = 7,   color = "grey30"),
      axis.ticks       = element_line(color = "grey40", linewidth = 0.20),
      axis.line        = element_line(color = "grey40", linewidth = 0.25),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid       = element_blank(),
      panel.border     = element_blank(),
      plot.background  = element_rect(fill = "white", color = NA),
      plot.margin      = margin(4, 6, 4, 4)
    )
}

qq_list <- list()

for (i in seq_along(hormones)) {
  qq_list[[paste0("PSQI~",horm_labs[i])]] <-
    make_qq(fits_psqi[[hormones[i]]]$resid,
            paste0("PSQI ~ ", horm_labs[i]),
            fits_psqi[[hormones[i]]]$type)
  qq_list[[paste0("PVS~",horm_labs[i])]] <-
    make_qq(fits_pvs[[hormones[i]]]$resid,
            paste0("PVS ~ ", horm_labs[i]),
            fits_pvs[[hormones[i]]]$type)
  qq_list[[paste0("MOD~",horm_labs[i])]] <-
    make_qq(fits_mod[[hormones[i]]]$resid,
            paste0("PVS ~ PSQI\u00d7", horm_labs[i]),
            fits_mod[[hormones[i]]]$type)
}

fig_qq <- cowplot::plot_grid(plotlist = qq_list, ncol = 4)

ggsave("qq_all_models.png", fig_qq,
       width = 14, height = 9, dpi = 300, bg = "white")

# ---- STEP 3: SLEEP-PVS COUPLING (delta + permutation) ----

# Delta variables
nat_delta <- nat %>%
  arrange(session) %>%
  mutate(
    d_pvs  = c(NA, diff(total_pvs_volume)),
    d_psqi = c(NA, diff(psqi))
  ) %>%
  filter(!is.na(d_pvs) & !is.na(d_psqi))

oc_delta_full <- oc %>%
  arrange(session) %>%
  mutate(
    d_pvs  = c(NA, diff(total_pvs_volume)),
    d_psqi = c(NA, diff(psqi))
  ) %>%
  filter(!is.na(d_pvs) & !is.na(d_psqi))

# Delta analyses also exclude the day following the outlier
oc_delta <- oc_delta_full %>%
  filter(!(session %in% outlier_session) &
           !(session %in% (outlier_session + 1)))

# Circular moving-block permutation of x
block_permute <- function(x, block_size) {
  n <- length(x)
  x_ext <- c(x, x)  # circular extension to avoid edge truncation
  starts <- sample(seq_len(n), size = ceiling(n / block_size), replace = TRUE)
  blocks <- lapply(starts, function(s) x_ext[s:(s + block_size - 1)])
  unlist(blocks)[seq_len(n)]
}

# Spearman delta correlation with block permutation
run_delta_cor <- function(delta_df, condition_lab, pvs_col = "d_pvs",
                          block_size = DEFAULT_BLOCK_SIZE, n_perm = 10000,
                          label = "PVS") {
  set.seed(42)  # reproducibility lock (see run_condition_perm)
  x <- delta_df$d_psqi
  y <- delta_df[[pvs_col]]
  
  ct    <- cor.test(x, y, method = "spearman", exact = FALSE)
  r_obs <- ct$estimate
  
  r_perm <- replicate(n_perm,
                      cor(block_permute(x, block_size), y,
                          method = "spearman", use = "complete.obs"))
  p_perm <- (sum(abs(r_perm) >= abs(r_obs)) + 1) / (n_perm + 1)  # +1/+1 continuity correction
  
  list(r = r_obs, p_spearman = ct$p.value,
       p_perm = p_perm, condition = condition_lab, label = label,
       block_size = block_size)
}

res_delta_nat <- run_delta_cor(nat_delta, "Physiological menstrual cycle",
                               label = "Total PVS")
res_delta_oc  <- run_delta_cor(oc_delta,  "OCP (primary, outlier day excluded)",
                               label = "Total PVS")
res_delta_nat
res_delta_oc

# Region-specific delta correlation (CSO/BG, physiological cycle only)
cso_col <- "pvs_volume_centrumsemiovale"
bg_col  <- "pvs_volume_basalganglia"

nat_delta <- nat_delta %>%
  mutate(
    d_pvs_cso = c(NA, diff(nat[[cso_col]]))[match(session, nat$session)],
    d_pvs_bg  = c(NA, diff(nat[[bg_col]]))[match(session, nat$session)]
  )
oc_delta <- oc_delta %>%
  mutate(
    d_pvs_cso = c(NA, diff(oc[[cso_col]]))[match(session, oc$session)]
  )

res_delta_nat_cso <- run_delta_cor(nat_delta, "Physiological menstrual cycle",
                                   pvs_col = "d_pvs_cso", label = "CSO")

res_delta_nat_bg <- run_delta_cor(nat_delta, "Physiological menstrual cycle",
                                  pvs_col = "d_pvs_bg", label = "BG")
res_delta_nat_cso
res_delta_nat_bg

# r_OC_CSO: not reported standalone, used below only
ct_oc_cso <- cor.test(oc_delta$d_psqi, oc_delta$d_pvs_cso,
                      method = "spearman", exact = FALSE)
res_delta_oc_cso <- list(r = ct_oc_cso$estimate, p_spearman = ct_oc_cso$p.value)

# Collinearity dCSO/dBG (physiological cycle only)
cor_cso_bg_nat <- cor.test(nat_delta$d_pvs_cso, nat_delta$d_pvs_bg,
                           method = "spearman", exact = FALSE)
cor_cso_bg_nat

# Partial Spearman correlations (physiological cycle only), permutation-based
partial_spearman <- function(x, y, z) {
  rx <- rank(x); ry <- rank(y); rz <- rank(z)
  resx <- residuals(lm(rx ~ rz))
  resy <- residuals(lm(ry ~ rz))
  cor(resx, resy, method = "pearson")  # Pearson on rank residuals = partial Spearman
}

run_partial_cor <- function(delta_df, condition_lab, x_col, y_col, z_col,
                            block_size = DEFAULT_BLOCK_SIZE, n_perm = 10000,
                            label = "") {
  set.seed(42)  # reproducibility lock (see run_condition_perm)
  x <- delta_df[[x_col]]; y <- delta_df[[y_col]]; z <- delta_df[[z_col]]
  r_obs <- partial_spearman(x, y, z)
  
  r_perm <- replicate(n_perm, partial_spearman(block_permute(x, block_size), y, z))
  p_perm <- (sum(abs(r_perm) >= abs(r_obs)) + 1) / (n_perm + 1)
  
  list(r = r_obs, p_perm = p_perm, condition = condition_lab, label = label)
}

res_partial_nat_cso <- run_partial_cor(nat_delta, "Physiological menstrual cycle",
                                       "d_psqi", "d_pvs_cso", "d_pvs_bg",
                                       label = "dPSQI-dCSO | dBG")
res_partial_nat_bg  <- run_partial_cor(nat_delta, "Physiological menstrual cycle",
                                       "d_psqi", "d_pvs_bg", "d_pvs_cso",
                                       label = "dPSQI-dBG | dCSO")
res_partial_nat_cso
res_partial_nat_bg

# Formal test of H0: r_physiol = r_OCP
test_coupling_difference <- function(delta_df1, delta_df2,
                                     cond1_lab = "physiol.", cond2_lab = "OCP",
                                     pvs_col = "d_pvs", label = "PVS",
                                     block_size = DIFF_TEST_BLOCK_SIZE,
                                     n_perm = 10000) {
  set.seed(42)  # reproducibility lock (see run_condition_perm)
  r1_obs <- cor(delta_df1$d_psqi, delta_df1[[pvs_col]],
                method = "spearman", use = "complete.obs")
  r2_obs <- cor(delta_df2$d_psqi, delta_df2[[pvs_col]],
                method = "spearman", use = "complete.obs")
  diff_obs <- r1_obs - r2_obs
  
  n1 <- nrow(delta_df1); n2 <- nrow(delta_df2)
  combined <- rbind(
    data.frame(d_psqi = delta_df1$d_psqi, d_pvs = delta_df1[[pvs_col]]),
    data.frame(d_psqi = delta_df2$d_psqi, d_pvs = delta_df2[[pvs_col]])
  )
  n_total  <- nrow(combined)
  block_id <- ceiling(seq_len(n_total) / block_size)
  n_blocks <- max(block_id)
  
  diff_perm <- replicate(n_perm, {
    block_order   <- sample(seq_len(n_blocks))
    reordered_idx <- unlist(lapply(block_order, function(b) which(block_id == b)))
    reordered     <- combined[reordered_idx, ]
    
    g1 <- reordered[1:n1, ]
    g2 <- reordered[(n1 + 1):n_total, ]
    r1 <- cor(g1$d_psqi, g1$d_pvs, method = "spearman", use = "complete.obs")
    r2 <- cor(g2$d_psqi, g2$d_pvs, method = "spearman", use = "complete.obs")
    r1 - r2
  })
  
  p_perm <- (sum(abs(diff_perm) >= abs(diff_obs)) + 1) / (n_perm + 1)
  
  list(r1 = r1_obs, r2 = r2_obs, diff = diff_obs, p_perm = p_perm, label = label)
}

diff_total <- test_coupling_difference(nat_delta, oc_delta, label = "Total PVS")

# Only one test in this family -> BH correction is the identity here.
p_diff_fdr <- p.adjust(diff_total$p_perm, method = "BH")
diff_total
p_diff_fdr

# Spread/dynamic-range summary for dPSQI (Table 1)
describe_spread <- function(x) {
  list(sd = sd(x, na.rm = TRUE), range = diff(range(x, na.rm = TRUE)),
       min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE),
       iqr = IQR(x, na.rm = TRUE), n = sum(!is.na(x)))
}
spread_stats <- list(
  physiol = describe_spread(nat_delta$d_psqi),
  ocp     = describe_spread(oc_delta$d_psqi)
)
spread_stats

# Variance-ratio test via block permutation (not parametric var.test())
run_variance_ratio_perm <- function(x1, x2, label,
                                    block_size = DIFF_TEST_BLOCK_SIZE,
                                    n_perm = 10000) {
  set.seed(42)  # reproducibility lock (see run_condition_perm)
  n1 <- length(x1); n2 <- length(x2)
  obs_ratio <- var(x1) / var(x2)
  
  combined <- c(x1, x2)
  n_total  <- length(combined)
  block_id <- ceiling(seq_len(n_total) / block_size)
  n_blocks <- max(block_id)
  
  ratio_perm <- replicate(n_perm, {
    block_order   <- sample(seq_len(n_blocks))
    reordered_idx <- unlist(lapply(block_order, function(b) which(block_id == b)))
    reordered     <- combined[reordered_idx]
    g1 <- reordered[1:n1]
    g2 <- reordered[(n1 + 1):n_total]
    var(g1) / var(g2)
  })
  
  # Log ratio for a symmetric two-sided test (ratios are right-skewed)
  obs_stat  <- abs(log(obs_ratio))
  perm_stat <- abs(log(ratio_perm))
  p_perm <- (sum(perm_stat >= obs_stat) + 1) / (n_perm + 1)
  
  list(ratio = obs_ratio, p_perm = p_perm)
}

vt_psqi <- run_variance_ratio_perm(nat_delta$d_psqi, oc_delta$d_psqi, "dPSQI (physiol. vs. OCP)")
vt_pvs  <- run_variance_ratio_perm(nat_delta$d_pvs,  oc_delta$d_pvs,  "dPVSV total (physiol. vs. OCP)")
vt_psqi
vt_pvs

# Mean comparison of the delta values themselves (not their correlation)
run_condition_perm_meandiff <- function(x1, x2, label,
                                        block_size = DIFF_TEST_BLOCK_SIZE,
                                        n_perm = 10000) {
  set.seed(42)  # reproducibility lock (see run_condition_perm)
  n1 <- length(x1); n2 <- length(x2)
  obs_diff <- mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)
  
  combined <- c(x1, x2)
  n_total  <- length(combined)
  block_id <- ceiling(seq_len(n_total) / block_size)
  n_blocks <- max(block_id)
  
  perm_diffs <- replicate(n_perm, {
    block_order   <- sample(seq_len(n_blocks))
    reordered_idx <- unlist(lapply(block_order, function(b) which(block_id == b)))
    x_reordered   <- combined[reordered_idx]
    mean(x_reordered[1:n1], na.rm = TRUE) -
      mean(x_reordered[(n1 + 1):n_total], na.rm = TRUE)
  })
  p_perm <- (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1)
  
  ci <- block_bootstrap_ci_meandiff(x1, x2,
                                    block_size1 = DEFAULT_BLOCK_SIZE,
                                    block_size2 = DEFAULT_BLOCK_SIZE)
  
  list(diff = obs_diff, ci_lo = ci["lo"], ci_hi = ci["hi"],
       p_perm = p_perm, label = label)
}

res_meandiff_dpvs  <- run_condition_perm_meandiff(nat_delta$d_pvs,  oc_delta$d_pvs,
                                                  "d PVSV, total")
res_meandiff_dpsqi <- run_condition_perm_meandiff(nat_delta$d_psqi, oc_delta$d_psqi,
                                                  "d PSQI total score")
res_meandiff_dpvs
res_meandiff_dpsqi

# Range-restriction sensitivity check (distribution-free)
oc_psqi_range <- range(oc_delta$d_psqi)

nat_restricted <- nat_delta %>%
  filter(d_psqi >= oc_psqi_range[1] & d_psqi <= oc_psqi_range[2])

run_range_restriction_boot <- function(nat_restricted, pvs_col, n_target,
                                       r_nat_full, r_oc_observed, label,
                                       n_boot = 5000) {
  set.seed(42)  # reproducibility lock (see run_condition_perm)
  if (nrow(nat_restricted) < 5) {
    return(invisible(NULL))
  }
  r_boot <- replicate(n_boot, {
    idx <- sample(seq_len(nrow(nat_restricted)), size = n_target, replace = TRUE)
    bd  <- nat_restricted[idx, ]
    if (sd(bd$d_psqi) == 0 || sd(bd[[pvs_col]]) == 0) return(NA)
    cor(bd$d_psqi, bd[[pvs_col]], method = "spearman")
  })
  r_boot <- r_boot[!is.na(r_boot)]
  
  dist_to_oc  <- abs(r_boot - r_oc_observed)
  dist_to_nat <- abs(r_boot - r_nat_full)
  pct_closer_to_oc <- mean(dist_to_oc < dist_to_nat) * 100
  
  invisible(list(r_boot = r_boot, pct_closer_to_oc = pct_closer_to_oc))
}

boot_total <- run_range_restriction_boot(nat_restricted, "d_pvs",
                                         nrow(oc_delta),
                                         res_delta_nat$r, res_delta_oc$r,
                                         label = "Total PVS")
boot_cso   <- run_range_restriction_boot(nat_restricted, "d_pvs_cso",
                                         nrow(oc_delta),
                                         res_delta_nat_cso$r, res_delta_oc_cso$r,
                                         label = "CSO")
boot_total
boot_cso

# 95% CIs via block bootstrap (10,000 iterations)
block_bootstrap_ci <- function(delta_df, pvs_col = "d_pvs",
                               block_size = DEFAULT_BLOCK_SIZE,
                               n_boot = 10000, conf = 0.95) {
  set.seed(42)  # reproducibility lock (see run_condition_perm)
  n <- nrow(delta_df)
  x <- delta_df$d_psqi
  y <- delta_df[[pvs_col]]
  block_id <- ceiling(seq_len(n) / block_size)
  n_blocks <- max(block_id)
  
  r_boot <- replicate(n_boot, {
    boot_blocks   <- sample(seq_len(n_blocks), size = n_blocks, replace = TRUE)
    boot_idx      <- unlist(lapply(boot_blocks, function(b) which(block_id == b)))
    boot_idx      <- boot_idx[seq_len(n)]  # truncate/pad to original length
    if (length(unique(x[boot_idx])) < 2 || length(unique(y[boot_idx])) < 2) return(NA)
    cor(x[boot_idx], y[boot_idx], method = "spearman", use = "complete.obs")
  })
  r_boot <- r_boot[!is.na(r_boot)]
  alpha <- 1 - conf
  c(lo = as.numeric(quantile(r_boot, alpha / 2)),
    hi = as.numeric(quantile(r_boot, 1 - alpha / 2)))
}

ci_nat      <- block_bootstrap_ci(nat_delta, "d_pvs",     DEFAULT_BLOCK_SIZE)
ci_oc       <- block_bootstrap_ci(oc_delta,  "d_pvs",     DEFAULT_BLOCK_SIZE)
ci_nat_cso  <- block_bootstrap_ci(nat_delta, "d_pvs_cso", DEFAULT_BLOCK_SIZE)
ci_nat_bg   <- block_bootstrap_ci(nat_delta, "d_pvs_bg",  DEFAULT_BLOCK_SIZE)
ci_nat
ci_oc
ci_nat_cso
ci_nat_bg

# Block-bootstrap CIs for the partial correlations
block_bootstrap_ci_partial <- function(delta_df, x_col, y_col, z_col,
                                       block_size = DEFAULT_BLOCK_SIZE,
                                       n_boot = 10000, conf = 0.95) {
  set.seed(42)  # reproducibility lock (see run_condition_perm)
  n <- nrow(delta_df)
  x <- delta_df[[x_col]]; y <- delta_df[[y_col]]; z <- delta_df[[z_col]]
  block_id <- ceiling(seq_len(n) / block_size)
  n_blocks <- max(block_id)
  
  r_boot <- replicate(n_boot, {
    boot_blocks <- sample(seq_len(n_blocks), size = n_blocks, replace = TRUE)
    boot_idx    <- unlist(lapply(boot_blocks, function(b) which(block_id == b)))
    boot_idx    <- boot_idx[seq_len(n)]
    xb <- x[boot_idx]; yb <- y[boot_idx]; zb <- z[boot_idx]
    if (length(unique(xb)) < 2 || length(unique(yb)) < 2 || length(unique(zb)) < 2) return(NA)
    tryCatch(partial_spearman(xb, yb, zb), error = function(e) NA)
  })
  r_boot <- r_boot[!is.na(r_boot)]
  alpha <- 1 - conf
  c(lo = as.numeric(quantile(r_boot, alpha / 2)),
    hi = as.numeric(quantile(r_boot, 1 - alpha / 2)))
}

ci_partial_cso <- block_bootstrap_ci_partial(nat_delta, "d_psqi", "d_pvs_cso", "d_pvs_bg")
ci_partial_bg  <- block_bootstrap_ci_partial(nat_delta, "d_psqi", "d_pvs_bg",  "d_pvs_cso")
ci_partial_cso
ci_partial_bg

# Lagged cross-correlations (lag +1 to +3)
run_lagged_cor <- function(delta_df, condition_lab, lag,
                           pvs_col = "d_pvs",
                           block_size = DEFAULT_BLOCK_SIZE, n_perm = 10000,
                           n_boot = 10000, conf = 0.95,
                           label = "PVS") {
  set.seed(42)  # reproducibility lock (see run_condition_perm)
  n <- nrow(delta_df)
  x <- delta_df$d_psqi[1:(n - lag)]                 # dPSQI_t
  y <- delta_df[[pvs_col]][(1 + lag):n]              # dPVS_{t+lag}
  n_pair <- length(x)
  
  ct    <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  r_obs <- ct$estimate
  
  r_perm <- replicate(n_perm,
                      cor(block_permute(x, block_size), y,
                          method = "spearman", use = "complete.obs"))
  p_perm <- (sum(abs(r_perm) >= abs(r_obs)) + 1) / (n_perm + 1)
  
  # Block-bootstrap CI, (x,y) pairs resampled jointly
  block_id_pair <- ceiling(seq_len(n_pair) / block_size)
  n_blocks_pair <- max(block_id_pair)
  r_boot <- replicate(n_boot, {
    boot_blocks <- sample(seq_len(n_blocks_pair), size = n_blocks_pair, replace = TRUE)
    boot_idx    <- unlist(lapply(boot_blocks, function(b) which(block_id_pair == b)))
    boot_idx    <- boot_idx[seq_len(n_pair)]
    if (length(unique(x[boot_idx])) < 2 || length(unique(y[boot_idx])) < 2) return(NA)
    cor(x[boot_idx], y[boot_idx], method = "spearman", use = "complete.obs")
  })
  r_boot <- r_boot[!is.na(r_boot)]
  alpha  <- 1 - conf
  ci_lo  <- as.numeric(quantile(r_boot, alpha / 2))
  ci_hi  <- as.numeric(quantile(r_boot, 1 - alpha / 2))
  
  list(r = r_obs, ci_lo = ci_lo, ci_hi = ci_hi, p_perm = p_perm, lag = lag,
       condition = condition_lab, label = label, n = n_pair)
}

lag_results <- list()
for (cond_df_name in c("nat_delta", "oc_delta")) {
  cond_df  <- get(cond_df_name)
  cond_lab <- ifelse(cond_df_name == "nat_delta",
                     "Physiological menstrual cycle", "OCP")
  for (lag in 1:3) {
    key <- paste0(cond_df_name, "_lag", lag)
    lag_results[[key]] <- run_lagged_cor(cond_df, cond_lab, lag,
                                         label = "Total PVS")
  }
}
# FDR across lag tests per condition (separate family from the lag-0 test)
p_lags_nat <- sapply(1:3, function(l) lag_results[[paste0("nat_delta_lag", l)]]$p_perm)
p_lags_oc  <- sapply(1:3, function(l) lag_results[[paste0("oc_delta_lag", l)]]$p_perm)
p_lags_nat_fdr <- p.adjust(p_lags_nat, "BH")
p_lags_oc_fdr  <- p.adjust(p_lags_oc,  "BH")
lag_results
p_lags_nat
p_lags_oc
p_lags_nat_fdr
p_lags_oc_fdr