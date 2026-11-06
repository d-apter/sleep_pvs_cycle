# ══════════════════════════════════════════════════════════════════════════════
# Statistical analysis
#   0. Data preparation (log-transformation of hormones)
#   1. Outlier exclusion (|z| > 3 for PVSV or 1-dPSQI)
#   2. Descriptive statistics
#   3. Condition comparison (physiological menstrual cycle vs. OCP)
#   4. Hormone regression and interaction models (physiological cycle)
#   5. Diagnostics (Durbin-Watson, Shapiro-Wilk, Q-Q plots, Cook's distance)
#   6. Day-to-day association of sleep quality and PVSV
#   7. Sensitivity analyses (GLS AR(1), GAM trend of study day, Johnson-Neyman,
#      rank-based GLS AR(1), block lengths 2-5)
# Requires a data frame `df` (see README) and optionally `psqi_scored`
# (1-dPSQI components). Results are printed as tables and the robustness
# results are written to CSV files in the working directory.
# ══════════════════════════════════════════════════════════════════════════════

library(dplyr)
library(nlme)
library(mgcv)
library(lmtest)
library(ggplot2)
library(cowplot)
filter  <- dplyr::filter
select  <- dplyr::select
arrange <- dplyr::arrange
mutate  <- dplyr::mutate

options(width = 160)
set.seed(42)

# Block length rule of thumb L ~ n^(1/3) (Hall, Horowitz & Jing, 1995):
# 3 days within one condition (and for bootstrap CIs), 4 days for tests on
# the combined dataset; sensitivity analyses with 2-5 days.
DEFAULT_BLOCK_SIZE   <- 3
DIFF_TEST_BLOCK_SIZE <- 4
BLOCK_SIZES_SENS     <- 2:5
GAM_K                <- 10
N_PERM               <- 10000
N_BOOT               <- 10000

# ══════════════════════════════════════════════════════════════════════════════
# Helper functions
# ══════════════════════════════════════════════════════════════════════════════
section <- function(title) cat(sprintf("\n==== %s ====\n", title))

# Prints a results table with numeric columns rounded to 3 significant digits
report <- function(x, title = NULL) {
  if (!is.null(title)) cat(sprintf("\n  %s\n", title))
  num <- vapply(x, is.numeric, logical(1))
  x[num] <- lapply(x[num], signif, 3)
  print(x, row.names = FALSE)
  invisible(x)
}

zs <- function(x) as.numeric(scale(x))

# Circular block resampling: x is split into blocks of consecutive values
# that are reassembled (circularly) with random starting points. Preserves
# the within-series autocorrelation of x but destroys its association with
# another series. Applied to seq_len(n), it returns block-bootstrap indices.
block_permute <- function(x, block_size) {
  n <- length(x)
  x_ext <- c(x, x)
  starts <- sample(seq_len(n), size = ceiling(n / block_size), replace = TRUE)
  blocks <- lapply(starts, function(s) x_ext[s:(s + block_size - 1)])
  unlist(blocks)[seq_len(n)]
}

block_boot_ci <- function(n, stat_fun, block_size = DEFAULT_BLOCK_SIZE,
                          n_boot = N_BOOT, conf = 0.95) {
  set.seed(42)
  stats <- replicate(n_boot, {
    idx <- block_permute(seq_len(n), block_size)
    tryCatch(stat_fun(idx), error = function(e) NA_real_)
  })
  alpha <- 1 - conf
  quantile(stats, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE)
}

spearman_safe <- function(x, y) {
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y, method = "spearman")
}

# ══════════════════════════════════════════════════════════════════════════════
# 0. Data preparation
# ══════════════════════════════════════════════════════════════════════════════
# Condition values in the data are "natural" and "OC"
add_log_hormones <- function(d) {
  d %>% mutate(
    log_E2  = log(`Estradiol (LCMS; pg/mL)`   + 1),
    log_P4  = log(`Progesterone (LCMS; ng/mL)` + 1),
    log_LH  = log(`LH (mIU/ml)`  + 1),
    log_FSH = log(`FSH (mIU/ml)` + 1)
  )
}

nat <- df %>%
  filter(condition == "natural") %>%
  arrange(session) %>%
  mutate(time = seq_len(n())) %>%
  add_log_hormones()

oc <- df %>%
  filter(condition == "OC") %>%
  arrange(session) %>%
  mutate(time = seq_len(n())) %>%
  add_log_hormones()

hormones  <- c("log_E2", "log_P4", "log_LH", "log_FSH")
horm_labs <- c("Estradiol", "Progesterone", "LH", "FSH")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Outlier exclusion (|z| > 3, predefined criterion)
# ══════════════════════════════════════════════════════════════════════════════
# Outlier days are excluded from all analyses. For the day-to-day analyses,
# both differences involving the outlier day are removed.
section("1. Outliers (|z| > 3)")

outlier_sessions <- function(d) {
  z <- abs(scale(d$total_pvs_volume)[, 1]) > 3 | abs(scale(d$psqi)[, 1]) > 3
  if (any(z)) report(d[z, c("condition", "session", "total_pvs_volume", "psqi")])
  d$session[z]
}

if (length(outlier_sessions(nat)) > 0) {
  stop("Outlier(s) in the physiological cycle: exclusion is implemented for the OCP condition only.")
}
outlier_session <- outlier_sessions(oc)
oc_clean <- oc %>% filter(!(session %in% outlier_session))

# ══════════════════════════════════════════════════════════════════════════════
# 2. Descriptive statistics
# ══════════════════════════════════════════════════════════════════════════════
section("2. Descriptive statistics")

desc_vars <- c("PVSV, total (mm3)"    = "total_pvs_volume",
               "PVSV, CSO (mm3)"      = "pvs_volume_centrumsemiovale",
               "PVSV, BG (mm3)"       = "pvs_volume_basalganglia",
               "1-dPSQI total score"  = "psqi",
               "Estradiol (pg/mL)"    = "Estradiol (LCMS; pg/mL)",
               "Progesterone (ng/mL)" = "Progesterone (LCMS; ng/mL)",
               "LH (mIU/mL)"          = "LH (mIU/ml)",
               "FSH (mIU/mL)"         = "FSH (mIU/ml)")

report(do.call(rbind, lapply(names(desc_vars), function(lab) {
  v <- desc_vars[[lab]]
  data.frame(variable = lab,
             physiol_mean = mean(nat[[v]]), physiol_sd = sd(nat[[v]]),
             ocp_mean = mean(oc_clean[[v]]), ocp_sd = sd(oc_clean[[v]]))
})), sprintf("Mean and SD (physiological cycle n = %d, OCP n = %d)", nrow(nat), nrow(oc_clean)))

# ══════════════════════════════════════════════════════════════════════════════
# 3. Condition comparison
# ══════════════════════════════════════════════════════════════════════════════
section("3. Condition comparison (physiological cycle vs. OCP)")

# Block permutation for two-sample statistics: both series are concatenated
# and split into fixed, non-overlapping blocks; the block order is shuffled
# and the first n1 values form the pseudo physiological-cycle group.
perm_two_sample <- function(x1, x2,
                            stat = function(a, b) mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
                            block_size = DIFF_TEST_BLOCK_SIZE, n_perm = N_PERM) {
  set.seed(42)
  n1 <- length(x1)
  x  <- c(x1, x2)
  n  <- length(x)
  block_id <- ceiling(seq_len(n) / block_size)
  n_blocks <- max(block_id)
  obs  <- stat(x1, x2)
  perm <- replicate(n_perm, {
    xr <- x[unlist(lapply(sample(seq_len(n_blocks)), function(b) which(block_id == b)))]
    stat(xr[seq_len(n1)], xr[(n1 + 1):n])
  })
  list(stat = obs, p_perm = (sum(abs(perm) >= abs(obs)) + 1) / (n_perm + 1))
}

# Block-bootstrap CI for the mean difference; each group is resampled separately
block_boot_ci_meandiff <- function(x1, x2, block_size = DEFAULT_BLOCK_SIZE,
                                   n_boot = N_BOOT, conf = 0.95) {
  set.seed(42)
  diffs <- replicate(n_boot, {
    i1 <- block_permute(seq_along(x1), block_size)
    i2 <- block_permute(seq_along(x2), block_size)
    mean(x1[i1], na.rm = TRUE) - mean(x2[i2], na.rm = TRUE)
  })
  alpha <- 1 - conf
  quantile(diffs, c(alpha / 2, 1 - alpha / 2), names = FALSE)
}

compare_conditions <- function(x1, x2, label, block_size = DIFF_TEST_BLOCK_SIZE,
                               ci = TRUE) {
  pt  <- perm_two_sample(x1, x2, block_size = block_size)
  ci_ <- if (ci) block_boot_ci_meandiff(x1, x2) else c(NA_real_, NA_real_)
  data.frame(variable = label, diff = pt$stat, ci_lo = ci_[1], ci_hi = ci_[2],
             p_perm = pt$p_perm, stringsAsFactors = FALSE)
}

cond_vars_hormone <- list(
  list(v = "Estradiol (LCMS; pg/mL)",    lab = "Estradiol"),
  list(v = "Progesterone (LCMS; ng/mL)", lab = "Progesterone"),
  list(v = "LH (mIU/ml)",                lab = "LH"),
  list(v = "FSH (mIU/ml)",               lab = "FSH")
)

cond_vars_outcome <- list(
  list(v = "total_pvs_volume",            lab = "PVS volume (total)"),
  list(v = "pvs_volume_centrumsemiovale", lab = "PVS volume (CSO)"),
  list(v = "pvs_volume_basalganglia",     lab = "PVS volume (BG)"),
  list(v = "psqi",                        lab = "PSQI global score")
)

run_family <- function(vars) {
  do.call(rbind, lapply(vars, function(v)
    compare_conditions(nat[[v$v]], oc_clean[[v$v]], v$lab)))
}

# Hormones are compared uncorrected (characterization of the hormonal
# profile); PVSV and 1-dPSQI are BH-corrected.
res_family1 <- run_family(cond_vars_hormone)
report(res_family1, "Hormones (difference physiological cycle - OCP, uncorrected)")

res_family2 <- run_family(cond_vars_outcome)
res_family2$p_fdr <- p.adjust(res_family2$p_perm, method = "BH")
report(res_family2, "PVSV and 1-dPSQI (BH-corrected)")

# 1-dPSQI components require `psqi_scored` (columns condition, session_num,
# C1_quality ... C7_daytime, optionally PSQI_total). C6 (sleep medication)
# is excluded, as it was zero on all days.
subcomp_vars <- list(
  list(v = "C1_quality",      lab = "C1 Subjective quality"),
  list(v = "C2_latency",      lab = "C2 Sleep latency"),
  list(v = "C3_duration",     lab = "C3 Sleep duration"),
  list(v = "C4_efficiency",   lab = "C4 Sleep efficiency"),
  list(v = "C5_disturbances", lab = "C5 Sleep disturbances"),
  list(v = "C7_daytime",      lab = "C7 Daytime dysfunction")
)

if (exists("psqi_scored")) {
  # session_num may use the session IDs of df or count days within condition
  ps_oc_ids <- psqi_scored$session_num[psqi_scored$condition == "OC"]
  if (all(ps_oc_ids %in% oc$session)) {
    out_ids <- outlier_session
  } else if (setequal(ps_oc_ids, oc$time)) {
    out_ids <- oc$time[oc$session %in% outlier_session]
  } else {
    stop("psqi_scored: session_num of the OCP condition matches neither the session IDs nor the days within condition.")
  }
  ps <- psqi_scored %>%
    filter(!(condition == "OC" & session_num %in% out_ids)) %>%
    arrange(factor(condition, levels = c("natural", "OC")), session_num)
  if (sum(ps$condition == "natural") != nrow(nat) || sum(ps$condition == "OC") != nrow(oc_clean)) {
    stop("psqi_scored: number of days does not match df after outlier exclusion.")
  }
  # Alignment check: a shift between days would produce widespread mismatches
  if ("PSQI_total" %in% names(ps)) {
    d_tot <- ps$PSQI_total - c(nat$psqi, oc_clean$psqi)
    if (sum(d_tot != 0, na.rm = TRUE) > 0.5 * sum(!is.na(d_tot))) {
      stop("psqi_scored: PSQI_total differs from psqi on most days; sessions are likely misaligned.")
    }
  }
  res_family3 <- do.call(rbind, lapply(subcomp_vars, function(v)
    compare_conditions(ps[[v$v]][ps$condition == "natural"],
                       ps[[v$v]][ps$condition == "OC"], v$lab, ci = FALSE)))
  report(res_family3, "1-dPSQI components (uncorrected)")
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Regression models (physiological cycle only)
#   All variables are z-standardized (standardized betas; centred before
#   forming the interaction term).
# ══════════════════════════════════════════════════════════════════════════════
section("4. Regression models (physiological cycle, standardized betas)")

nat_z <- nat %>%
  mutate(psqi             = zs(psqi),
         total_pvs_volume = zs(total_pvs_volume),
         across(all_of(hormones), zs))

# OLS, or GLS with AR(1) errors if the Durbin-Watson test indicates residual
# autocorrelation
fit_ar <- function(formula_str, data) {
  m0 <- lm(as.formula(formula_str), data = data)
  dw <- dwtest(m0, alternative = "two.sided")
  if (dw$p.value < 0.05) {
    m <- gls(as.formula(formula_str), data = data,
             correlation = corAR1(form = ~time), method = "ML")
    list(model = m, type = "GLS AR(1)", dw_p = dw$p.value, dw_stat = dw$statistic,
         resid = residuals(m, type = "normalized"))
  } else {
    list(model = m0, type = "OLS", dw_p = dw$p.value, dw_stat = dw$statistic,
         resid = residuals(m0))
  }
}

get_coef <- function(fit, term) {
  m <- fit$model
  if (inherits(m, "gls")) {
    s  <- summary(m)$tTable
    b  <- s[term, "Value"]
    se <- s[term, "Std.Error"]
    p  <- s[term, "p-value"]
    ci <- b + c(-1, 1) * qt(0.975, m$dims$N - m$dims$p) * se
  } else {
    s  <- summary(m)$coefficients
    b  <- s[term, "Estimate"]
    se <- s[term, "Std. Error"]
    p  <- s[term, "Pr(>|t|)"]
    ci <- confint(m)[term, ]
  }
  list(beta = b, se = se, p = p, ci_lo = ci[1], ci_hi = ci[2])
}

# Model 1: 1-dPSQI ~ hormone; Model 2: PVSV ~ hormone;
# Model 3: PVSV ~ 1-dPSQI x hormone
fams <- c("1a" = "Model 1", "1b" = "Model 2", "2" = "Model 3")

fml_for <- function(fam, h, extra = "") {
  base <- switch(fam,
                 "1a" = paste0("psqi ~ ", h),
                 "1b" = paste0("total_pvs_volume ~ ", h),
                 "2"  = paste0("total_pvs_volume ~ psqi * ", h))
  paste0(base, extra)
}
term_for <- function(fam, h) if (fam == "2") paste0("psqi:", h) else h
label_for <- function(fam, lab) {
  switch(fam,
         "1a" = paste0("PSQI ~ ", lab),
         "1b" = paste0("PVS ~ ", lab),
         "2"  = paste0("PVS ~ PSQI x ", lab))
}

fits     <- list()
coef_tab <- list()
for (fam in names(fams)) {
  fits[[fam]] <- list()
  for (i in seq_along(hormones)) {
    h   <- hormones[i]
    fit <- fit_ar(fml_for(fam, h), nat_z)
    fits[[fam]][[h]] <- fit
    cf  <- get_coef(fit, term_for(fam, h))
    coef_tab[[length(coef_tab) + 1]] <- data.frame(
      fam = fam, hormone = h, label = label_for(fam, horm_labs[i]), type = fit$type,
      beta = cf$beta, se = cf$se, ci_lo = cf$ci_lo, ci_hi = cf$ci_hi, p = cf$p,
      stringsAsFactors = FALSE)
  }
}
coef_tab <- do.call(rbind, coef_tab)
# BH correction across the four hormones within each model
coef_tab$p_fdr <- ave(coef_tab$p, coef_tab$fam, FUN = function(p) p.adjust(p, "BH"))
report(coef_tab[, c("label", "type", "beta", "se", "ci_lo", "ci_hi", "p", "p_fdr")])

# ══════════════════════════════════════════════════════════════════════════════
# 5. Diagnostics
# ══════════════════════════════════════════════════════════════════════════════
section("5. Diagnostics")

cooks_thr <- 4 / nrow(nat_z)
diag_tab <- do.call(rbind, lapply(seq_len(nrow(coef_tab)), function(k) {
  fit <- fits[[coef_tab$fam[k]]][[coef_tab$hormone[k]]]
  cd  <- if (fit$type == "OLS") cooks.distance(fit$model) else NA_real_
  data.frame(model = coef_tab$label[k], dw = as.numeric(fit$dw_stat), dw_p = fit$dw_p,
             sw_p = shapiro.test(fit$resid)$p.value, max_cooks_d = max(cd),
             session_max_d = if (fit$type == "OLS") nat_z$session[which.max(cd)] else NA,
             n_above_4n = sum(cd > cooks_thr))
}))
report(diag_tab, "Durbin-Watson, Shapiro-Wilk and Cook's distance (threshold 4/n)")

# Leave-one-out for FDR-significant models: the day with the largest Cook's
# distance (if above 4/n) is removed and the model re-estimated after
# re-standardizing within the reduced sample.
run_loo <- function(fam, h, label) {
  fit <- fits[[fam]][[h]]
  if (fit$type != "OLS") return(NULL)
  cd <- cooks.distance(fit$model)
  if (max(cd) <= cooks_thr) return(NULL)
  sess  <- nat$session[which.max(cd)]
  fml   <- fml_for(fam, h)
  d_loo <- nat %>% filter(session != sess)
  for (v in all.vars(as.formula(fml))) d_loo[[v]] <- zs(d_loo[[v]])
  cf <- get_coef(fit_ar(fml, d_loo), term_for(fam, h))
  data.frame(model = label, excluded_session = sess, cooks_d = max(cd), n = nrow(d_loo),
             beta = cf$beta, se = cf$se, ci_lo = cf$ci_lo, ci_hi = cf$ci_hi, p = cf$p)
}

loo_tab <- do.call(rbind, lapply(which(coef_tab$p_fdr < 0.05), function(k)
  run_loo(coef_tab$fam[k], coef_tab$hormone[k], coef_tab$label[k])))
if (!is.null(loo_tab)) report(loo_tab, "Leave-one-out (most influential day)")

make_qq <- function(res, label) {
  qq <- qqnorm(res, plot.it = FALSE)
  ggplot(data.frame(theoretical = qq$x, sample = qq$y), aes(theoretical, sample)) +
    geom_abline(slope = 1, intercept = 0, color = "grey60", linewidth = 0.5, linetype = "dashed") +
    geom_point(size = 1.8, shape = 21, fill = "#3A8C8C", color = "white", stroke = 0.3) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 2.5,
             label = sprintf("Shapiro-Wilk p = %.3f", shapiro.test(res)$p.value)) +
    labs(title = label, x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_classic(base_size = 8) +
    theme(plot.title = element_text(face = "bold", size = 7.5))
}

qq_list <- lapply(seq_len(nrow(coef_tab)), function(k)
  make_qq(fits[[coef_tab$fam[k]]][[coef_tab$hormone[k]]]$resid, coef_tab$label[k]))
ggsave("qq_all_models.png", cowplot::plot_grid(plotlist = qq_list, ncol = 4),
       width = 14, height = 9, dpi = 300, bg = "white")

# ══════════════════════════════════════════════════════════════════════════════
# 6. Day-to-day association of sleep quality and PVSV
# ══════════════════════════════════════════════════════════════════════════════
section("6. Day-to-day association (first-order differences)")

cso_col <- "pvs_volume_centrumsemiovale"
bg_col  <- "pvs_volume_basalganglia"

make_delta <- function(d) {
  d %>%
    arrange(session) %>%
    mutate(d_pvs     = c(NA, diff(total_pvs_volume)),
           d_psqi    = c(NA, diff(psqi)),
           d_pvs_cso = c(NA, diff(.data[[cso_col]])),
           d_pvs_bg  = c(NA, diff(.data[[bg_col]]))) %>%
    filter(!is.na(d_pvs) & !is.na(d_psqi))
}

nat_delta <- make_delta(nat)
oc_delta  <- make_delta(oc) %>%
  filter(!(session %in% outlier_session) & !(session %in% (outlier_session + 1)))

# Autocorrelation of level and difference series (ACF lags 1-3, Ljung-Box)
autocorr_row <- function(x, label, max_lag = 3) {
  x  <- x[!is.na(x)]
  ac <- acf(x, lag.max = max_lag, plot = FALSE)$acf[-1]
  lb <- Box.test(x, lag = max_lag, type = "Ljung-Box")
  data.frame(series = label, r1 = ac[1], r2 = ac[2], r3 = ac[3],
             Q = as.numeric(lb$statistic), df = max_lag, p = lb$p.value)
}

acf_table <- rbind(
  autocorr_row(nat$psqi,                  "Physiol.: PSQI (level)"),
  autocorr_row(nat$total_pvs_volume,      "Physiol.: PVS (level)"),
  autocorr_row(nat_delta$d_psqi,          "Physiol.: dPSQI"),
  autocorr_row(nat_delta$d_pvs,           "Physiol.: dPVS"),
  autocorr_row(nat_delta$d_pvs_cso,       "Physiol.: dPVS CSO"),
  autocorr_row(nat_delta$d_pvs_bg,        "Physiol.: dPVS BG"),
  autocorr_row(oc_clean$psqi,             "OCP: PSQI (level)"),
  autocorr_row(oc_clean$total_pvs_volume, "OCP: PVS (level)"),
  autocorr_row(oc_delta$d_psqi,           "OCP: dPSQI"),
  autocorr_row(oc_delta$d_pvs,            "OCP: dPVS")
)
report(acf_table, "Autocorrelation")
write.csv(acf_table, "autocorrelation_series.csv", row.names = FALSE)

# Spearman correlation of the differences; p from circular block permutation
# of x (y held fixed), 95% CI from circular block bootstrap
coupling_test <- function(x, y, block_size = DEFAULT_BLOCK_SIZE, n_perm = N_PERM,
                          ci = TRUE) {
  set.seed(42)
  r_obs  <- cor(x, y, method = "spearman")
  r_perm <- replicate(n_perm, cor(block_permute(x, block_size), y, method = "spearman"))
  p_perm <- (sum(abs(r_perm) >= abs(r_obs)) + 1) / (n_perm + 1)
  ci_ <- if (ci) block_boot_ci(length(x), function(i) spearman_safe(x[i], y[i]), block_size)
         else c(NA_real_, NA_real_)
  list(r = r_obs, p_perm = p_perm, ci_lo = ci_[1], ci_hi = ci_[2])
}

# Partial Spearman correlation (rank-residual approach, equivalent to
# ppcor::pcor.test(..., method = "spearman")); x is block-permuted, z fixed
partial_spearman <- function(x, y, z) {
  rz <- rank(z)
  cor(residuals(lm(rank(x) ~ rz)), residuals(lm(rank(y) ~ rz)))
}

partial_test <- function(x, y, z, block_size = DEFAULT_BLOCK_SIZE, n_perm = N_PERM,
                         ci = TRUE) {
  set.seed(42)
  r_obs  <- partial_spearman(x, y, z)
  r_perm <- replicate(n_perm, partial_spearman(block_permute(x, block_size), y, z))
  p_perm <- (sum(abs(r_perm) >= abs(r_obs)) + 1) / (n_perm + 1)
  ci_ <- if (ci) block_boot_ci(length(x), function(i) {
    if (sd(x[i]) == 0 || sd(y[i]) == 0 || sd(z[i]) == 0) return(NA_real_)
    partial_spearman(x[i], y[i], z[i])
  }, block_size) else c(NA_real_, NA_real_)
  list(r = r_obs, p_perm = p_perm, ci_lo = ci_[1], ci_hi = ci_[2])
}

as_row <- function(res, label, n) {
  data.frame(test = label, n = n, rho = res$r, ci_lo = res$ci_lo, ci_hi = res$ci_hi,
             p_perm = res$p_perm)
}

res_delta_nat <- coupling_test(nat_delta$d_psqi, nat_delta$d_pvs)
res_delta_oc  <- coupling_test(oc_delta$d_psqi,  oc_delta$d_pvs)
delta_tab <- rbind(as_row(res_delta_nat, "Physiol. cycle, total PVS", nrow(nat_delta)),
                   as_row(res_delta_oc,  "OCP, total PVS",            nrow(oc_delta)))

# Regional and partial correlations for each condition with a significant
# total-PVSV association
delta_sets <- list(
  list(d = nat_delta, lab = "Physiol. cycle", res = res_delta_nat),
  list(d = oc_delta,  lab = "OCP",            res = res_delta_oc)
)
delta_sets_sig <- Filter(function(ds) ds$res$p_perm < 0.05, delta_sets)

regional <- list()
for (ds in delta_sets_sig) {
  d <- ds$d
  regional[[ds$lab]] <- list(
    cso         = coupling_test(d$d_psqi, d$d_pvs_cso),
    bg          = coupling_test(d$d_psqi, d$d_pvs_bg),
    partial_cso = partial_test(d$d_psqi, d$d_pvs_cso, d$d_pvs_bg),
    partial_bg  = partial_test(d$d_psqi, d$d_pvs_bg,  d$d_pvs_cso)
  )
  r <- regional[[ds$lab]]
  delta_tab <- rbind(delta_tab,
                     as_row(r$cso,         paste0(ds$lab, ", CSO"), nrow(d)),
                     as_row(r$bg,          paste0(ds$lab, ", BG"), nrow(d)),
                     as_row(r$partial_cso, paste0(ds$lab, ", partial dCSO | dBG"), nrow(d)),
                     as_row(r$partial_bg,  paste0(ds$lab, ", partial dBG | dCSO"), nrow(d)))
}

# Difference between conditions (Gelman & Stern, 2006): (dPSQI, dPVS) pairs
# stay together; blocks of consecutive days are reordered and split into two
# pseudo-groups of the original sizes.
test_coupling_difference <- function(d1, d2, pvs_col = "d_pvs",
                                     block_size = DIFF_TEST_BLOCK_SIZE, n_perm = N_PERM) {
  set.seed(42)
  r1 <- cor(d1$d_psqi, d1[[pvs_col]], method = "spearman")
  r2 <- cor(d2$d_psqi, d2[[pvs_col]], method = "spearman")
  comb <- rbind(data.frame(x = d1$d_psqi, y = d1[[pvs_col]]),
                data.frame(x = d2$d_psqi, y = d2[[pvs_col]]))
  n1 <- nrow(d1)
  n  <- nrow(comb)
  block_id <- ceiling(seq_len(n) / block_size)
  n_blocks <- max(block_id)
  diff_perm <- replicate(n_perm, {
    rc <- comb[unlist(lapply(sample(seq_len(n_blocks)), function(b) which(block_id == b))), ]
    cor(rc$x[1:n1], rc$y[1:n1], method = "spearman") -
      cor(rc$x[(n1 + 1):n], rc$y[(n1 + 1):n], method = "spearman")
  })
  list(r1 = r1, r2 = r2, diff = r1 - r2,
       p_perm = (sum(abs(diff_perm) >= abs(r1 - r2)) + 1) / (n_perm + 1))
}

diff_total <- test_coupling_difference(nat_delta, oc_delta)
delta_tab <- rbind(delta_tab, data.frame(
  test = "Difference physiol. - OCP (delta rho)", n = NA, rho = diff_total$diff,
  ci_lo = NA, ci_hi = NA, p_perm = diff_total$p_perm))
report(delta_tab, "Spearman correlations of day-to-day changes")

# Variance ratio of the day-to-day changes (block permutation of the log
# ratio, as var.test() assumes independent observations)
log_var_ratio <- function(a, b) log(var(a) / var(b))
vr_vars <- c(dPSQI = "d_psqi", `dPVS total` = "d_pvs")
report(do.call(rbind, lapply(names(vr_vars), function(lab) {
  col <- vr_vars[[lab]]
  vr  <- perm_two_sample(nat_delta[[col]], oc_delta[[col]], stat = log_var_ratio)
  data.frame(variable = lab, sd_physiol = sd(nat_delta[[col]]), sd_ocp = sd(oc_delta[[col]]),
             variance_ratio = exp(vr$stat), p_perm = vr$p_perm)
})), "Variance ratio of day-to-day changes (physiological cycle / OCP)")

# Mean day-to-day change by condition
res_family4 <- rbind(
  compare_conditions(nat_delta$d_pvs,  oc_delta$d_pvs,  "dPVS total"),
  compare_conditions(nat_delta$d_psqi, oc_delta$d_psqi, "dPSQI total score")
)
report(res_family4, "Mean day-to-day change (uncorrected)")

# Range restriction (Sackett & Yang, 2000): physiological-cycle differences
# within the dPSQI range observed under OCP, bootstrapped to the OCP sample size
oc_psqi_range  <- range(oc_delta$d_psqi)
nat_restricted <- nat_delta %>%
  filter(d_psqi >= oc_psqi_range[1] & d_psqi <= oc_psqi_range[2])

if (nrow(nat_restricted) >= 5) {
  set.seed(42)
  r_boot <- replicate(5000, {
    bd <- nat_restricted[sample(seq_len(nrow(nat_restricted)), nrow(oc_delta), replace = TRUE), ]
    spearman_safe(bd$d_psqi, bd$d_pvs)
  })
  r_boot <- r_boot[!is.na(r_boot)]
  report(data.frame(n_within_range = nrow(nat_restricted), median_rho = median(r_boot),
                    pi_lo = quantile(r_boot, 0.025, names = FALSE),
                    pi_hi = quantile(r_boot, 0.975, names = FALSE),
                    rho_physiol = res_delta_nat$r, rho_ocp = res_delta_oc$r),
         sprintf("Range restriction (OCP dPSQI range %g to %g, 5,000 bootstrap replicates)",
                 oc_psqi_range[1], oc_psqi_range[2]))
}

# Lagged cross-correlations: dPSQI_t with dPVS_{t+lag} (directional
# assumption that sleep quality affects PVSV). Pairs are matched by session,
# so that no pair spans the days excluded around the OCP outlier.
lag_pairs <- function(d, lag) {
  idx  <- match(d$session + lag, d$session)
  keep <- !is.na(idx)
  list(x = d$d_psqi[keep], y = d$d_pvs[idx[keep]])
}

lag_results <- do.call(rbind, lapply(list(list(d = nat_delta, lab = "Physiol. cycle"),
                                          list(d = oc_delta,  lab = "OCP")), function(cond) {
  res <- do.call(rbind, lapply(1:3, function(lag) {
    lp <- lag_pairs(cond$d, lag)
    r  <- coupling_test(lp$x, lp$y)
    data.frame(condition = cond$lab, lag = lag, n = length(lp$x), rho = r$r,
               ci_lo = r$ci_lo, ci_hi = r$ci_hi, p_perm = r$p_perm)
  }))
  res$p_fdr <- p.adjust(res$p_perm, method = "BH")
  res
}))
report(lag_results, "Lagged cross-correlations (BH across lags within condition)")

# ══════════════════════════════════════════════════════════════════════════════
# 7. Sensitivity analyses
# ══════════════════════════════════════════════════════════════════════════════
section("7. Sensitivity analyses")

fit_gls_ar1 <- function(fml, data) {
  list(model = gls(as.formula(fml), data = data,
                   correlation = corAR1(form = ~time), method = "ML"),
       type = "GLS AR(1)")
}

corr_param <- function(m) {
  cs <- m$modelStruct$corStruct
  if (is.null(cs)) NA_real_ else as.numeric(coef(cs, unconstrained = FALSE))
}

# Regression models refitted as GLS with AR(1) errors regardless of the
# Durbin-Watson results
rob_reg <- do.call(rbind, lapply(seq_len(nrow(coef_tab)), function(k) {
  fam <- coef_tab$fam[k]
  h   <- coef_tab$hormone[k]
  g   <- fit_gls_ar1(fml_for(fam, h), nat_z)
  cg  <- get_coef(g, term_for(fam, h))
  data.frame(model = coef_tab$label[k],
             beta = coef_tab$beta[k], ci_lo = coef_tab$ci_lo[k], ci_hi = coef_tab$ci_hi[k],
             p = coef_tab$p[k], p_fdr = coef_tab$p_fdr[k],
             gls_ar1_beta = cg$beta, gls_ar1_ci_lo = cg$ci_lo, gls_ar1_ci_hi = cg$ci_hi,
             gls_ar1_p = cg$p, gls_ar1_phi = corr_param(g$model),
             stringsAsFactors = FALSE)
}))

# Smooth trend of study day (GAM, all fitted by ML so that AIC is comparable):
#   hormone: primary model
#   trend:   hormone term(s) replaced by a penalized smooth of study day
#   both:    primary model plus the smooth of study day
trend_fml_for <- function(fam) {
  sm <- sprintf("s(time, k = %d)", GAM_K)
  switch(fam,
         "1a" = paste0("psqi ~ ", sm),
         "1b" = paste0("total_pvs_volume ~ ", sm),
         "2"  = paste0("total_pvs_volume ~ psqi + ", sm))
}

gam_coef <- function(m, term) {
  pt <- summary(m)$p.table
  b  <- pt[term, "Estimate"]
  se <- pt[term, "Std. Error"]
  ci <- b + c(-1, 1) * qt(0.975, df.residual(m)) * se
  list(beta = b, se = se, p = pt[term, "Pr(>|t|)"], ci_lo = ci[1], ci_hi = ci[2])
}

gam_res <- do.call(rbind, lapply(seq_len(nrow(coef_tab)), function(k) {
  fam  <- coef_tab$fam[k]
  h    <- coef_tab$hormone[k]
  m_h  <- mgcv::gam(as.formula(fml_for(fam, h)), data = nat_z, method = "ML")
  m_t  <- mgcv::gam(as.formula(trend_fml_for(fam)), data = nat_z, method = "ML")
  m_b  <- mgcv::gam(as.formula(fml_for(fam, h, sprintf(" + s(time, k = %d)", GAM_K))),
                    data = nat_z, method = "ML")
  cb   <- gam_coef(m_b, term_for(fam, h))
  s_t  <- summary(m_t)$s.table
  s_b  <- summary(m_b)$s.table
  data.frame(model = coef_tab$label[k],
             aic_hormone = AIC(m_h), aic_trend = AIC(m_t), aic_both = AIC(m_b),
             trend_edf = s_t[1, "edf"], trend_p = s_t[1, "p-value"],
             both_beta = cb$beta, both_ci_lo = cb$ci_lo, both_ci_hi = cb$ci_hi, both_p = cb$p,
             both_trend_edf = s_b[1, "edf"], both_trend_p = s_b[1, "p-value"],
             stringsAsFactors = FALSE)
}))
rob_reg <- merge(rob_reg, gam_res, by = "model", sort = FALSE)
report(rob_reg[, c("model", "beta", "p_fdr", "gls_ar1_beta", "gls_ar1_p", "gls_ar1_phi",
                   "aic_hormone", "aic_trend", "both_beta", "both_p")],
       "Regression models: GLS AR(1) and smooth trend of study day")
write.csv(rob_reg, "robustness_regression.csv", row.names = FALSE)

# Simple slopes (-1 SD, mean, +1 SD) and Johnson-Neyman regions for the
# gonadotropin interactions (t-based critical value; raw hormone units)
johnson_neyman <- function(m, h, data, log_raw, label, n_grid = 1000) {
  b  <- coef(m)
  V  <- vcov(m)
  t1 <- "psqi"
  t3 <- paste0("psqi:", h)
  tc <- qt(0.975, df.residual(m))
  slope_at <- function(z) as.numeric(b[t1] + b[t3] * z)
  se_at    <- function(z) as.numeric(sqrt(V[t1, t1] + z^2 * V[t3, t3] + 2 * z * V[t1, t3]))
  to_raw   <- function(z) exp(z * sd(log_raw) + mean(log_raw)) - 1

  z <- c(-1, 0, 1)
  slopes <- data.frame(hormone = label, at = c("-1 SD", "mean", "+1 SD"), level = to_raw(z),
                       slope = slope_at(z), ci_lo = slope_at(z) - tc * se_at(z),
                       ci_hi = slope_at(z) + tc * se_at(z),
                       p = 2 * pt(-abs(slope_at(z) / se_at(z)), df.residual(m)))

  grid <- seq(min(data[[h]]), max(data[[h]]), length.out = n_grid)
  tval <- slope_at(grid) / se_at(grid)
  runs <- rle(abs(tval) > tc)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1
  sig <- which(runs$values)
  t_days <- slope_at(data[[h]]) / se_at(data[[h]])
  regions <- data.frame(hormone = label,
                        direction = ifelse(tval[starts[sig]] > 0, "positive", "negative"),
                        from = to_raw(grid[starts[sig]]), to = to_raw(grid[ends[sig]]),
                        days_sig_pos = sum(t_days > tc), days_sig_neg = sum(t_days < -tc))
  list(slopes = slopes, regions = regions)
}

jn <- lapply(c("log_LH", "log_FSH"), function(h)
  johnson_neyman(fits[["2"]][[h]]$model, h, nat_z, nat[[h]], horm_labs[hormones == h]))
report(do.call(rbind, lapply(jn, `[[`, "slopes")), "Simple slopes of PVSV on 1-dPSQI")
report(do.call(rbind, lapply(jn, `[[`, "regions")), "Johnson-Neyman regions of significance")

# Day-to-day association: rank-based GLS with AR(1) errors
prep_rank_delta <- function(d, cond) {
  data.frame(condition = cond, time = d$time,
             x = zs(rank(d$d_psqi)), y = zs(rank(d$d_pvs)))
}
d_rank <- rbind(prep_rank_delta(nat_delta, "natural"), prep_rank_delta(oc_delta, "OC"))
d_rank$condition <- factor(d_rank$condition, levels = c("natural", "OC"))

fit_ar1_delta <- function(data, fml, grouped = FALSE) {
  cs <- if (grouped) corAR1(form = ~ time | condition) else corAR1(form = ~ time)
  list(model = gls(as.formula(fml), data = data, correlation = cs, method = "ML"),
       type = "GLS AR(1)")
}

ar1_nat <- fit_ar1_delta(subset(d_rank, condition == "natural"), "y ~ x")
ar1_oc  <- fit_ar1_delta(subset(d_rank, condition == "OC"),      "y ~ x")
ar1_int <- fit_ar1_delta(d_rank, "y ~ x * condition", grouped = TRUE)
c_nat <- get_coef(ar1_nat, "x")
c_oc  <- get_coef(ar1_oc,  "x")
c_int <- get_coef(ar1_int, "x:conditionOC")

rob_coupling <- data.frame(
  effect   = c("dPSQI-dPVS, physiol.", "dPSQI-dPVS, OCP", "Difference (physiol. - OCP)"),
  estimate = c(res_delta_nat$r, res_delta_oc$r, diff_total$diff),
  ci_lo    = c(res_delta_nat$ci_lo, res_delta_oc$ci_lo, NA),
  ci_hi    = c(res_delta_nat$ci_hi, res_delta_oc$ci_hi, NA),
  p_perm   = c(res_delta_nat$p_perm, res_delta_oc$p_perm, diff_total$p_perm),
  ar1_beta  = c(c_nat$beta,  c_oc$beta,  -c_int$beta),
  ar1_ci_lo = c(c_nat$ci_lo, c_oc$ci_lo, -c_int$ci_hi),
  ar1_ci_hi = c(c_nat$ci_hi, c_oc$ci_hi, -c_int$ci_lo),
  ar1_p     = c(c_nat$p,     c_oc$p,     c_int$p),
  ar1_phi   = c(corr_param(ar1_nat$model), corr_param(ar1_oc$model), corr_param(ar1_int$model)),
  stringsAsFactors = FALSE)
report(rob_coupling, "Day-to-day association: Spearman (block permutation) vs. rank-based GLS AR(1)")
write.csv(rob_coupling, "robustness_day_to_day.csv", row.names = FALSE)

# Block length: every block-permutation test repeated with 2-5 days, with
# the same correction as in the primary analysis (BH for PVSV/1-dPSQI and
# across lags; all other tests uncorrected)
run_block <- function(b) {
  rows <- list()
  add <- function(family, test, primary, p, p_fdr = NA_real_) {
    rows[[length(rows) + 1]] <<- data.frame(family = family, test = test,
                                            primary_block = primary, block = b,
                                            p = p, p_fdr = p_fdr, stringsAsFactors = FALSE)
  }
  cond_p <- function(x1, x2) compare_conditions(x1, x2, "", block_size = b, ci = FALSE)$p_perm

  for (v in cond_vars_hormone) {
    add("Family 1: hormones", v$lab, DIFF_TEST_BLOCK_SIZE, cond_p(nat[[v$v]], oc_clean[[v$v]]))
  }
  p2 <- sapply(cond_vars_outcome, function(v) cond_p(nat[[v$v]], oc_clean[[v$v]]))
  f2 <- p.adjust(p2, method = "BH")
  for (k in seq_along(cond_vars_outcome)) {
    add("Family 2: PVS volume and PSQI", cond_vars_outcome[[k]]$lab, DIFF_TEST_BLOCK_SIZE, p2[k], f2[k])
  }
  if (exists("ps")) {
    for (v in subcomp_vars) {
      add("Family 3: PSQI sub-components", v$lab, DIFF_TEST_BLOCK_SIZE,
          cond_p(ps[[v$v]][ps$condition == "natural"], ps[[v$v]][ps$condition == "OC"]))
    }
  }
  add("Family 4: delta means", "dPVS total",        DIFF_TEST_BLOCK_SIZE, cond_p(nat_delta$d_pvs,  oc_delta$d_pvs))
  add("Family 4: delta means", "dPSQI total score", DIFF_TEST_BLOCK_SIZE, cond_p(nat_delta$d_psqi, oc_delta$d_psqi))
  for (v in list(list(col = "d_psqi", lab = "dPSQI"), list(col = "d_pvs", lab = "dPVS total"))) {
    add("Variance ratio", v$lab, DIFF_TEST_BLOCK_SIZE,
        perm_two_sample(nat_delta[[v$col]], oc_delta[[v$col]], stat = log_var_ratio, block_size = b)$p_perm)
  }

  cp <- function(x, y) coupling_test(x, y, block_size = b, ci = FALSE)$p_perm
  add("Delta correlations", "Physiol. cycle, total PVS", DEFAULT_BLOCK_SIZE, cp(nat_delta$d_psqi, nat_delta$d_pvs))
  add("Delta correlations", "OCP, total PVS",            DEFAULT_BLOCK_SIZE, cp(oc_delta$d_psqi,  oc_delta$d_pvs))
  for (ds in delta_sets_sig) {
    d <- ds$d
    add("Delta correlations", paste0(ds$lab, ", CSO"), DEFAULT_BLOCK_SIZE, cp(d$d_psqi, d$d_pvs_cso))
    add("Delta correlations", paste0(ds$lab, ", BG"),  DEFAULT_BLOCK_SIZE, cp(d$d_psqi, d$d_pvs_bg))
    add("Delta correlations", paste0(ds$lab, ", partial dCSO | dBG"), DEFAULT_BLOCK_SIZE,
        partial_test(d$d_psqi, d$d_pvs_cso, d$d_pvs_bg, block_size = b, ci = FALSE)$p_perm)
    add("Delta correlations", paste0(ds$lab, ", partial dBG | dCSO"), DEFAULT_BLOCK_SIZE,
        partial_test(d$d_psqi, d$d_pvs_bg, d$d_pvs_cso, block_size = b, ci = FALSE)$p_perm)
  }
  add("Coupling difference", "Physiol. - OCP, total PVS", DIFF_TEST_BLOCK_SIZE,
      test_coupling_difference(nat_delta, oc_delta, block_size = b)$p_perm)

  for (cond in list(list(d = nat_delta, lab = "Physiol. cycle"), list(d = oc_delta, lab = "OCP"))) {
    pl <- sapply(1:3, function(lag) {
      lp <- lag_pairs(cond$d, lag)
      cp(lp$x, lp$y)
    })
    fl <- p.adjust(pl, method = "BH")
    for (lag in 1:3) {
      add("Lagged correlations", sprintf("%s, lag +%d", cond$lab, lag), DEFAULT_BLOCK_SIZE, pl[lag], fl[lag])
    }
  }
  do.call(rbind, rows)
}

blk_all <- do.call(rbind, lapply(BLOCK_SIZES_SENS, run_block))
blk_key <- factor(paste(blk_all$family, blk_all$test), levels = unique(paste(blk_all$family, blk_all$test)))

rob_block <- do.call(rbind, lapply(split(blk_all, blk_key), function(d) {
  prim <- d[d$block == d$primary_block[1], ]
  fdr  <- d$p_fdr[!is.na(d$p_fdr)]
  data.frame(family = d$family[1], test = d$test[1], primary_block = d$primary_block[1],
             p_primary = prim$p, p_min = min(d$p), p_max = max(d$p),
             p_fdr_primary = prim$p_fdr,
             p_fdr_min = if (length(fdr)) min(fdr) else NA_real_,
             p_fdr_max = if (length(fdr)) max(fdr) else NA_real_,
             stringsAsFactors = FALSE)
}))
rownames(rob_block) <- NULL
report(rob_block[, -1], sprintf("Block-length sensitivity (%d-%d days)",
                                min(BLOCK_SIZES_SENS), max(BLOCK_SIZES_SENS)))
write.csv(rob_block, "robustness_block_length.csv", row.names = FALSE)
