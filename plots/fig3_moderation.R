# Figure 3: gonadotropin trajectories, conditional slopes and
# Johnson-Neyman interaction plots
#   a  FSH and LH across the physiological menstrual cycle
#   b  Conditional slope of PVSV on 1-dPSQI across the cycle (FSH model)
#   c  Conditional slope of PVSV on 1-dPSQI across the cycle (LH model)
#   d  PVSV ~ 1-dPSQI x FSH (Johnson-Neyman)
#   e  PVSV ~ 1-dPSQI x LH (Johnson-Neyman)
# Requires `df`; run from the repository root or from plots/.

source(if (file.exists("fig_style.R")) "fig_style.R" else file.path("plots", "fig_style.R"))

FIG_WIDTH  <- 180   # mm
FIG_HEIGHT <- 185   # mm
N_JN_LINES <- 41    # regression lines per Johnson-Neyman panel
jn_gradient <- list(FSH = c("#9CCFC0", "#1F4D42"), LH = c("#E6C48F", "#5E4116"))

nat <- prep_nat(df)
nat <- nat[order(nat$cycle_day), ]

# ---- MODELS ----
fit_int <- function(h) lm(as.formula(paste0("total_pvs_volume ~ psqi * ", h)), data = nat)
fits <- lapply(c(log_E2 = "log_E2", log_P4 = "log_P4", log_LH = "log_LH", log_FSH = "log_FSH"), fit_int)
p_int <- sapply(names(fits), function(h)
  summary(fits[[h]])$coefficients[paste0("psqi:", h), "Pr(>|t|)"])
p_int_fdr <- p.adjust(p_int, method = "BH")

std_beta <- function(h) {
  d <- nat %>% mutate(across(c(total_pvs_volume, psqi, all_of(h)), ~ as.numeric(scale(.x))))
  coef(lm(as.formula(paste0("total_pvs_volume ~ psqi * ", h)), data = d))[paste0("psqi:", h)]
}
stat_label <- function(h) {
  sprintf("beta == %.2f*','~~italic(p)[FDR]~'%s'", std_beta(h), fmt_p(p_int_fdr[h]))
}

# Conditional slope of PVSV on 1-dPSQI at moderator value m (t-based CI)
cond_slope <- function(fit, h, m) {
  b <- coef(fit); V <- vcov(fit); it <- paste0("psqi:", h)
  tc <- qt(1 - ALPHA_SIG / 2, df.residual(fit))
  s  <- b["psqi"] + b[it] * m
  se <- sqrt(V["psqi", "psqi"] + m^2 * V[it, it] + 2 * m * V["psqi", it])
  data.frame(slope = as.numeric(s), lo = as.numeric(s - tc * se),
             hi = as.numeric(s + tc * se), sig = as.numeric(abs(s / se)) > tc)
}

x_cycle <- scale_x_continuous(breaks = c(1, seq(5, 30, 5)), limits = c(1, 30),
                              expand = expansion(add = 0.3))

# ---- PANEL A: TRAJECTORIES ----
smooth_series <- function(x, y, n = 300) {
  s <- spline(x = x, y = y, n = n)
  data.frame(x = s$x, y = s$y)
}
traj <- bind_rows(
  smooth_series(nat$cycle_day, as.numeric(scale(nat$log_FSH))) %>% mutate(hormone = "FSH"),
  smooth_series(nat$cycle_day, as.numeric(scale(nat$log_LH)))  %>% mutate(hormone = "LH")
)

panel_a <- ggplot(traj, aes(x = x, y = y)) +
  band_layers(phase_bands) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = LW_THIN, linetype = "22") +
  geom_ribbon(aes(ymin = 0, ymax = y, fill = hormone), alpha = 0.18, show.legend = FALSE) +
  geom_line(aes(color = hormone, linetype = hormone), linewidth = LW_MAIN) +
  scale_color_manual(values = hormone_colors, name = NULL) +
  scale_fill_manual(values = hormone_colors) +
  scale_linetype_manual(values = hormone_linetype, name = NULL) +
  x_cycle +
  scale_y_continuous(breaks = seq(-2, 4, 2)) +
  coord_cartesian(ylim = c(-2.2, 5.8), clip = "off") +
  labs(x = "Cycle day", y = "z-score (log)") +
  theme_fig() +
  theme(legend.position = c(0.99, 0.62), legend.justification = c(1, 0.5),
        legend.key.width = unit(18, "pt"))

# ---- PANELS B/C: CONDITIONAL SLOPES ACROSS THE CYCLE ----
slope_panel <- function(h, label) {
  cs <- cbind(nat[, c("cycle_day", "phase")], cond_slope(fits[[h]], h, nat[[h]]))
  sm <- data.frame(x  = spline(cs$cycle_day, cs$slope, n = 300)$x,
                   y  = spline(cs$cycle_day, cs$slope, n = 300)$y,
                   lo = spline(cs$cycle_day, cs$lo,    n = 300)$y,
                   hi = spline(cs$cycle_day, cs$hi,    n = 300)$y)
  col <- hormone_colors[label]
  y_top <- max(cs$hi) + 0.12 * diff(range(c(cs$lo, cs$hi)))
  ggplot() +
    band_layers(phase_bands) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = LW_THIN, linetype = "22") +
    geom_ribbon(data = sm, aes(x = x, ymin = lo, ymax = hi), fill = col, alpha = 0.2) +
    geom_line(data = sm, aes(x = x, y = lo), color = col, linewidth = LW_THIN) +
    geom_line(data = sm, aes(x = x, y = hi), color = col, linewidth = LW_THIN) +
    geom_line(data = sm, aes(x = x, y = y), color = col, linewidth = LW_MAIN) +
    geom_point(data = cs, aes(x = cycle_day, y = slope, fill = phase, shape = phase),
               size = POINT_SIZE, stroke = POINT_STROKE, color = "black") +
    tag_layer(paste(label, "model")) +
    scale_fill_manual(values = phase_colors, drop = FALSE) +
    scale_shape_manual(values = phase_shapes, drop = FALSE) +
    x_cycle +
    coord_cartesian(ylim = c(min(cs$lo), y_top), clip = "off") +
    labs(x = "Cycle day", y = "Simple slope") +
    theme_fig() +
    theme(legend.position = "none")
}
panel_b <- slope_panel("log_FSH", "FSH")
panel_c <- slope_panel("log_LH",  "LH")

# ---- PANELS D/E: JOHNSON-NEYMAN INTERACTION PLOTS ----
# One regression line of PVSV on 1-dPSQI per moderator value: coloured and
# solid where the conditional slope is significant (t-based), grey and
# dotted where it is not.
jn_panel <- function(h, label) {
  fit   <- fits[[h]]
  b     <- coef(fit)
  m_seq <- seq(min(nat[[h]]), max(nat[[h]]), length.out = N_JN_LINES)
  cs    <- cond_slope(fit, h, m_seq)
  x_rng <- range(nat$psqi)
  lines_df <- do.call(rbind, lapply(seq_along(m_seq), function(k) data.frame(
    id = k, x = x_rng, m_raw = exp(m_seq[k]) - 1,
    slope_type = ifelse(cs$sig[k], "p < 0.05", "n.s."),
    y = b["(Intercept)"] + b[h] * m_seq[k] + cs$slope[k] * x_rng)))
  lines_df$slope_type <- factor(lines_df$slope_type, levels = c("p < 0.05", "n.s."))

  y_rng <- range(nat$total_pvs_volume)
  y_lim <- y_rng + c(-0.12, 0.22) * diff(y_rng)

  ggplot() +
    geom_line(data = subset(lines_df, slope_type == "n.s."),
              aes(x = x, y = y, group = id, linetype = slope_type),
              color = "grey55", linewidth = LW_THIN) +
    geom_line(data = subset(lines_df, slope_type == "p < 0.05"),
              aes(x = x, y = y, group = id, color = m_raw, linetype = slope_type),
              linewidth = LW_THIN) +
    geom_point(data = nat, aes(x = psqi, y = total_pvs_volume, fill = phase, shape = phase),
               size = POINT_SIZE, stroke = POINT_STROKE, color = "black") +
    stat_layer(stat_label(h), x = x_rng[1], y = y_lim[2]) +
    scale_color_gradient(low = jn_gradient[[label]][1], high = jn_gradient[[label]][2],
                         trans = "log", breaks = breaks_log(n = 4),
                         limits = range(lines_df$m_raw),
                         labels = label_number(accuracy = 1),
                         name = paste0(label, " (mIU/mL)"),
                         guide = guide_colorbar(order = 1, barheight = unit(40, "pt"),
                                                barwidth = unit(6, "pt"))) +
    scale_linetype_manual(values = c("p < 0.05" = "solid", "n.s." = "11"),
                          name = "Slope", drop = FALSE,
                          guide = guide_legend(order = 2, override.aes = list(
                            color = c("black", "grey55"), linewidth = 0.4))) +
    scale_fill_manual(values = phase_colors, guide = "none") +
    scale_shape_manual(values = phase_shapes, guide = "none") +
    scale_x_continuous(breaks = pretty_breaks(n = 5)) +
    coord_cartesian(ylim = y_lim) +
    labs(x = "1-dPSQI", y = expression("PVSV (mm"^3*")")) +
    theme_fig() +
    theme(legend.position = "right", legend.key.width = unit(14, "pt"),
          legend.spacing.y = unit(4, "pt"), legend.box.spacing = unit(2, "pt"))
}
panel_d <- jn_panel("log_FSH", "FSH")
panel_e <- jn_panel("log_LH",  "LH")

# ---- ASSEMBLE ----
left_col <- label_panels(panel_a, panel_b, panel_c, ncol = 1, labels = c("a", "b", "c"),
                         align = "v", axis = "lr")
row_de   <- label_panels(panel_d, panel_e, ncol = 2, labels = c("d", "e"),
                         align = "h", axis = "tb")
fig3 <- cowplot::plot_grid(left_col, row_de,
                           phase_legend(phase_colors, phase_shapes, "Cycle phase"),
                           ncol = 1, rel_heights = c(1.45, 1, 0.07))
save_fig(fig3, "Fig3", FIG_WIDTH, FIG_HEIGHT)

# Johnson-Neyman regions for the figure legend (raw hormone units)
for (h in c("log_FSH", "log_LH")) {
  m_grid <- seq(min(nat[[h]]), max(nat[[h]]), length.out = 2000)
  cs <- cond_slope(fits[[h]], h, m_grid)
  r  <- rle(cs$sig); e <- cumsum(r$lengths); s <- e - r$lengths + 1
  for (k in which(r$values)) {
    cat(sprintf("%s: significant %s slope for %.2f to %.2f mIU/mL\n", h,
                ifelse(cs$slope[s[k]] > 0, "positive", "negative"),
                exp(m_grid[s[k]]) - 1, exp(m_grid[e[k]]) - 1))
  }
}
