# Figure 2: hormone associations with sleep quality and PVSV
#   a  1-dPSQI ~ LH
#   b  PVSV ~ estradiol
# Figure 4: daily dynamics of sleep quality and PVSV
#   a  Trajectories, physiological menstrual cycle
#   b  Trajectories, OCP
#   c  Day-to-day changes, physiological menstrual cycle
#   d  Day-to-day changes, OCP
# Statistics mirror the main analysis script (OLS with z-standardized
# variables and BH correction within model groups; delta Spearman with
# circular block permutation, block length 3, 10,000 permutations, seed 42).
# Requires `df`; run from the repository root or from plots/.

source(if (file.exists("fig_style.R")) "fig_style.R" else file.path("plots", "fig_style.R"))

N_PERM     <- 10000
BLOCK_SIZE <- 3

nat <- prep_nat(df)
oc  <- prep_oc(df)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2
# ══════════════════════════════════════════════════════════════════════════════
hormones <- c("log_E2", "log_P4", "log_LH", "log_FSH")
zs <- function(x) as.numeric(scale(x))
nat_z <- nat %>% mutate(psqi = zs(psqi), total_pvs_volume = zs(total_pvs_volume),
                        across(all_of(hormones), zs))

reg_stats <- function(outcome) {
  res <- t(sapply(hormones, function(h) {
    s <- summary(lm(as.formula(paste(outcome, "~", h)), data = nat_z))$coefficients
    c(beta = s[h, "Estimate"], p = s[h, "Pr(>|t|)"])
  }))
  data.frame(hormone = hormones, beta = res[, "beta"], p = res[, "p"],
             p_fdr = p.adjust(res[, "p"], method = "BH"))
}
st_psqi <- reg_stats("psqi")
st_pvs  <- reg_stats("total_pvs_volume")

reg_label <- function(st, h) {
  r <- st[st$hormone == h, ]
  sprintf("beta == %.2f*','~~italic(p)[FDR]~'%s'", r$beta, fmt_p(r$p_fdr))
}

reg_panel <- function(x_var, y_var, x_lab, y_lab, color, label) {
  fit <- lm(as.formula(paste(y_var, "~", x_var)), data = nat)
  x_seq <- seq(min(nat[[x_var]]), max(nat[[x_var]]), length.out = 200)
  pr <- predict(fit, newdata = setNames(data.frame(x_seq), x_var), interval = "confidence")
  reg_df <- data.frame(x = x_seq, fit = pr[, "fit"], lo = pr[, "lwr"], hi = pr[, "upr"])
  y_rng <- range(c(nat[[y_var]], reg_df$lo, reg_df$hi))
  y_top <- y_rng[2] + 0.18 * diff(y_rng)
  ggplot() +
    geom_ribbon(data = reg_df, aes(x = x, ymin = lo, ymax = hi), fill = color, alpha = 0.2) +
    geom_line(data = reg_df, aes(x = x, y = lo), color = color, linewidth = LW_THIN) +
    geom_line(data = reg_df, aes(x = x, y = hi), color = color, linewidth = LW_THIN) +
    geom_line(data = reg_df, aes(x = x, y = fit), color = color, linewidth = LW_MAIN) +
    geom_point(data = nat, aes(x = .data[[x_var]], y = .data[[y_var]], fill = phase, shape = phase),
               size = POINT_SIZE, stroke = POINT_STROKE, color = "black") +
    stat_layer(label, x = min(nat[[x_var]]), y = y_top) +
    scale_fill_manual(values = phase_colors, guide = "none") +
    scale_shape_manual(values = phase_shapes, guide = "none") +
    scale_x_continuous(breaks = pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = pretty_breaks(n = 5)) +
    coord_cartesian(ylim = c(y_rng[1], y_top)) +
    labs(x = x_lab, y = y_lab) +
    theme_fig()
}

f2_a <- reg_panel("log_LH", "psqi", "LH (log-transformed)", "1-dPSQI",
                  hormone_colors["LH"], reg_label(st_psqi, "log_LH"))
f2_b <- reg_panel("log_E2", "total_pvs_volume", "Estradiol (log-transformed)",
                  expression("PVSV (mm"^3*")"),
                  hormone_colors["Estradiol"], reg_label(st_pvs, "log_E2"))

fig2 <- cowplot::plot_grid(
  label_panels(f2_a, f2_b, ncol = 2, labels = c("a", "b")),
  phase_legend(phase_colors, phase_shapes, "Cycle phase"),
  ncol = 1, rel_heights = c(1, 0.1))
save_fig(fig2, "Fig2", 180, 80)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4
# ══════════════════════════════════════════════════════════════════════════════

# ---- Trajectories (z-scored across both conditions, visualization only) ----
traj_vars <- c(pvs = "total_pvs_volume", psqi = "psqi", Estradiol = "log_E2",
               Progesterone = "log_P4", LH = "log_LH", FSH = "log_FSH")
pooled <- lapply(traj_vars, function(v) {
  x <- c(nat[[v]], oc[[v]])
  list(m = mean(x, na.rm = TRUE), s = sd(x, na.rm = TRUE))
})
z_pooled <- function(d) {
  for (k in names(traj_vars)) d[[paste0(k, "_z")]] <- (d[[traj_vars[k]]] - pooled[[k]]$m) / pooled[[k]]$s
  d[order(d$cycle_day), ]
}
nat_p <- z_pooled(nat)
oc_p  <- z_pooled(oc)

loess_curve <- function(d, var, span = 0.35) {
  fit <- loess(as.formula(paste(var, "~ cycle_day")), data = d, span = span)
  x <- seq(min(d$cycle_day), max(d$cycle_day), length.out = 200)
  data.frame(x = x, y = predict(fit, newdata = data.frame(cycle_day = x)))
}

traj_ylim <- range(unlist(lapply(list(nat_p, oc_p), function(d)
  d[, paste0(c("Estradiol", "Progesterone", "LH", "FSH"), "_z")])), na.rm = TRUE)
traj_ylim <- c(floor(traj_ylim[1]), ceiling(traj_ylim[2]) + 1)

traj_panel <- function(d, bands) {
  long <- do.call(rbind, lapply(c("Estradiol", "Progesterone", "LH", "FSH"), function(h)
    data.frame(cycle_day = d$cycle_day, z = d[[paste0(h, "_z")]], hormone = h)))
  long$hormone <- factor(long$hormone, levels = c("Estradiol", "Progesterone", "LH", "FSH"))
  pvs_c  <- loess_curve(d, "pvs_z")
  psqi_c <- loess_curve(d, "psqi_z")

  # Inline labels at the right end, separated if they would overlap
  end_labels <- c("PVSV", "1-dPSQI")
  y_end <- c(tail(pvs_c$y, 1), tail(psqi_c$y, 1))
  gap <- 0.07 * diff(traj_ylim)
  if (abs(diff(y_end)) < gap) {
    y_end <- mean(y_end) + c(1, -1) * ifelse(y_end[1] >= y_end[2], 1, -1) * gap / 2
  }

  ggplot() +
    band_layers(bands) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = LW_THIN, linetype = "22") +
    geom_ribbon(data = long, aes(x = cycle_day, ymin = 0, ymax = z, fill = hormone,
                                 color = hormone, linetype = hormone),
                alpha = 0.25, linewidth = 0.35, outline.type = "upper") +
    geom_line(data = psqi_c, aes(x = x, y = y), color = psqi_color, linewidth = 0.8,
              linetype = "52") +
    geom_line(data = pvs_c, aes(x = x, y = y), color = pvs_color, linewidth = 0.8) +
    annotate("text", x = 30.6, y = y_end, label = end_labels, hjust = 0,
             size = pt_to_size(TXT), fontface = "bold", family = FONT_FAMILY,
             color = c(pvs_color, psqi_color)) +
    scale_fill_manual(values = hormone_colors, name = "Hormones") +
    scale_color_manual(values = hormone_colors, name = "Hormones") +
    scale_linetype_manual(values = traj_linetype, name = "Hormones") +
    scale_x_continuous(breaks = c(1, seq(5, 30, 5))) +
    coord_cartesian(xlim = c(0.7, 30.3), ylim = traj_ylim, expand = FALSE, clip = "off") +
    labs(x = "Cycle day", y = "z-score") +
    theme_fig() +
    theme(legend.position = "none", plot.margin = margin(12, 34, 4, 4))
}

f4_a <- traj_panel(nat_p, phase_bands)
f4_b <- traj_panel(oc_p,  phase_bands_oc)

leg_hormones <- extract_legend(
  f4_a + theme(legend.position = "bottom", legend.direction = "horizontal",
               legend.key.width = unit(20, "pt")))

# ---- Day-to-day changes ----
make_delta <- function(d) {
  d <- d[order(d$session), ]
  data.frame(session = d$session[-1], phase = d$phase[-1],
             d_pvs = diff(d$total_pvs_volume), d_psqi = diff(d$psqi))
}
outlier_oc <- oc$session[abs(scale(oc$psqi)[, 1]) > 3 | abs(scale(oc$total_pvs_volume)[, 1]) > 3]
nat_delta <- make_delta(nat)
oc_delta  <- make_delta(oc)
oc_delta  <- oc_delta[!(oc_delta$session %in% c(outlier_oc, outlier_oc + 1)), ]

block_permute <- function(x, block_size) {
  n <- length(x)
  x_ext <- c(x, x)
  starts <- sample(seq_len(n), size = ceiling(n / block_size), replace = TRUE)
  unlist(lapply(starts, function(s) x_ext[s:(s + block_size - 1)]))[seq_len(n)]
}
coupling <- function(d) {
  set.seed(42)
  r <- cor(d$d_psqi, d$d_pvs, method = "spearman")
  r_perm <- replicate(N_PERM, cor(block_permute(d$d_psqi, BLOCK_SIZE), d$d_pvs, method = "spearman"))
  list(r = r, p = (sum(abs(r_perm) >= abs(r)) + 1) / (N_PERM + 1))
}
cp_nat <- coupling(nat_delta)
cp_oc  <- coupling(oc_delta)
cp_label <- function(cp) sprintf("italic(rho) == %.2f*','~~italic(p)[perm]~'%s'", cp$r, fmt_p(cp$p))

x_lim_d <- range(c(nat_delta$d_psqi, oc_delta$d_psqi)) + c(-0.5, 0.5)
y_rng_d <- range(c(nat_delta$d_pvs, oc_delta$d_pvs))
y_lim_d <- c(y_rng_d[1] - 0.05 * diff(y_rng_d), y_rng_d[2] + 0.25 * diff(y_rng_d))

delta_panel <- function(d, cp, colors, shapes) {
  fit <- lm(d_pvs ~ d_psqi, data = d)
  x_seq <- seq(min(d$d_psqi), max(d$d_psqi), length.out = 200)
  pr <- predict(fit, newdata = data.frame(d_psqi = x_seq), interval = "confidence")
  reg_df <- data.frame(x = x_seq, fit = pr[, "fit"], lo = pr[, "lwr"], hi = pr[, "upr"])
  col <- "#4A8585"
  ggplot() +
    geom_hline(yintercept = 0, color = "grey40", linewidth = LW_THIN, linetype = "22") +
    geom_vline(xintercept = 0, color = "grey40", linewidth = LW_THIN, linetype = "22") +
    geom_ribbon(data = reg_df, aes(x = x, ymin = lo, ymax = hi), fill = col, alpha = 0.2) +
    geom_line(data = reg_df, aes(x = x, y = lo), color = col, linewidth = LW_THIN) +
    geom_line(data = reg_df, aes(x = x, y = hi), color = col, linewidth = LW_THIN) +
    geom_line(data = reg_df, aes(x = x, y = fit), color = col, linewidth = LW_MAIN) +
    geom_point(data = d, aes(x = d_psqi, y = d_pvs, fill = phase, shape = phase),
               size = POINT_SIZE, stroke = POINT_STROKE, color = "black") +
    stat_layer(cp_label(cp), x = x_lim_d[1], y = y_lim_d[2]) +
    scale_fill_manual(values = colors, guide = "none") +
    scale_shape_manual(values = shapes, guide = "none") +
    scale_x_continuous(breaks = pretty_breaks(n = 6)) +
    coord_cartesian(xlim = x_lim_d, ylim = y_lim_d, expand = FALSE) +
    labs(x = expression(Delta*"1-dPSQI"), y = expression(Delta*"PVSV (mm"^3*")")) +
    theme_fig()
}

f4_c <- delta_panel(nat_delta, cp_nat, phase_colors,    phase_shapes)
f4_d <- delta_panel(oc_delta,  cp_oc,  phase_colors_oc, phase_shapes_oc)

# Condition tags above the panels, right-aligned to the edge of each column,
# so that they line up between a/c and b/d and never touch the plotted data
condition_tags <- function(row) {
  tag <- function(label, x) grid::textGrob(
    label, x = unit(x, "npc") - unit(4, "pt"), y = unit(1, "npc") - unit(2, "pt"),
    hjust = 1, vjust = 1, gp = grid::gpar(fontsize = TXT, fontfamily = FONT_FAMILY))
  cowplot::ggdraw(row) +
    cowplot::draw_grob(tag("Physiological menstrual cycle", 0.5)) +
    cowplot::draw_grob(tag("OCP", 1))
}
fig4 <- cowplot::plot_grid(
  condition_tags(label_panels(f4_a, f4_b, ncol = 2, labels = c("a", "b"))),
  leg_hormones,
  condition_tags(label_panels(f4_c, f4_d, ncol = 2, labels = c("c", "d"))),
  cowplot::plot_grid(phase_legend(phase_colors, phase_shapes, "Cycle phase"),
                     phase_legend(phase_colors_oc, phase_shapes_oc, "OCP phase"),
                     ncol = 2),
  ncol = 1, rel_heights = c(1, 0.09, 1.15, 0.09))
save_fig(fig4, "Fig4", 180, 150)

cat(sprintf("Fig 4 day-to-day association: physiol. rho = %.3f, p_perm = %.4f; OCP rho = %.3f, p_perm = %.4f\n",
            cp_nat$r, cp_nat$p, cp_oc$r, cp_oc$p))
