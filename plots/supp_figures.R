# Supplementary figures 2-5 (Supplementary Fig. 1, the segmentation example,
# is made from the imaging data)
#   SuppFig2  Remaining hormone regressions
#             a-c 1-dPSQI ~ estradiol / progesterone / FSH
#             d-f PVSV ~ progesterone / LH / FSH
#   SuppFig3  Spearman correlations of each hormone with the 1-dPSQI
#             components, physiological cycle (a-d: E2, P4, LH, FSH)
#   SuppFig4  Spearman correlations of PVSV with the 1-dPSQI components
#             (a physiological cycle, b OCP)
#   SuppFig5  Day-to-day association by PVSV region (physiological cycle):
#             a CSO, b BG
# Statistics mirror the main analysis script. Requires `df` and `psqi_scored`;
# run from the repository root or from plots/.

source(if (file.exists("fig_style.R")) "fig_style.R" else file.path("plots", "fig_style.R"))

N_PERM     <- 10000
BLOCK_SIZE <- 3
# FALSE: PVSV-component correlations on daily levels (Supplementary Fig. 4);
# TRUE: on day-to-day changes
PVS_SPIDER_DELTA <- FALSE

nat <- prep_nat(df)
oc  <- prep_oc(df)
outlier_oc <- oc$session[abs(scale(oc$psqi)[, 1]) > 3 | abs(scale(oc$total_pvs_volume)[, 1]) > 3]
oc_clean <- oc[!(oc$session %in% outlier_oc), ]

# ══════════════════════════════════════════════════════════════════════════════
# Non-significant hormone regressions
# ══════════════════════════════════════════════════════════════════════════════
hormones <- c("log_E2", "log_P4", "log_LH", "log_FSH")
horm_lab <- c(log_E2 = "Estradiol", log_P4 = "Progesterone", log_LH = "LH", log_FSH = "FSH")
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

reg_panel <- function(x_var, y_var, y_lab, st) {
  color <- hormone_colors[horm_lab[x_var]]
  fit <- lm(as.formula(paste(y_var, "~", x_var)), data = nat)
  x_seq <- seq(min(nat[[x_var]]), max(nat[[x_var]]), length.out = 200)
  pr <- predict(fit, newdata = setNames(data.frame(x_seq), x_var), interval = "confidence")
  reg_df <- data.frame(x = x_seq, fit = pr[, "fit"], lo = pr[, "lwr"], hi = pr[, "upr"])
  y_rng <- range(c(nat[[y_var]], reg_df$lo, reg_df$hi))
  y_top <- y_rng[2] + 0.2 * diff(y_rng)
  ggplot() +
    geom_ribbon(data = reg_df, aes(x = x, ymin = lo, ymax = hi), fill = color, alpha = 0.2) +
    geom_line(data = reg_df, aes(x = x, y = lo), color = color, linewidth = LW_THIN) +
    geom_line(data = reg_df, aes(x = x, y = hi), color = color, linewidth = LW_THIN) +
    geom_line(data = reg_df, aes(x = x, y = fit), color = color, linewidth = LW_MAIN) +
    geom_point(data = nat, aes(x = .data[[x_var]], y = .data[[y_var]], fill = phase, shape = phase),
               size = POINT_SIZE, stroke = POINT_STROKE, color = "black") +
    stat_layer(reg_label(st, x_var), x = min(nat[[x_var]]), y = y_top) +
    scale_fill_manual(values = phase_colors, guide = "none") +
    scale_shape_manual(values = phase_shapes, guide = "none") +
    scale_x_continuous(breaks = pretty_breaks(n = 4)) +
    scale_y_continuous(breaks = pretty_breaks(n = 5)) +
    coord_cartesian(ylim = c(y_rng[1], y_top)) +
    labs(x = paste0(horm_lab[x_var], " (log-transformed)"), y = y_lab) +
    theme_fig()
}

pvs_lab <- expression("PVSV (mm"^3*")")
s_reg <- cowplot::plot_grid(
  label_panels(reg_panel("log_E2",  "psqi", "1-dPSQI", st_psqi),
               reg_panel("log_P4",  "psqi", "1-dPSQI", st_psqi),
               reg_panel("log_FSH", "psqi", "1-dPSQI", st_psqi),
               reg_panel("log_P4",  "total_pvs_volume", pvs_lab, st_pvs),
               reg_panel("log_LH",  "total_pvs_volume", pvs_lab, st_pvs),
               reg_panel("log_FSH", "total_pvs_volume", pvs_lab, st_pvs),
               ncol = 3, labels = letters[1:6]),
  phase_legend(phase_colors, phase_shapes, "Cycle phase"),
  ncol = 1, rel_heights = c(1, 0.06))
save_fig(s_reg, "SuppFig2", 180, 125)

# ══════════════════════════════════════════════════════════════════════════════
# Regional day-to-day association (physiological cycle)
# ══════════════════════════════════════════════════════════════════════════════
make_delta <- function(d) {
  d <- d[order(d$session), ]
  data.frame(session = d$session[-1], phase = d$phase[-1],
             d_psqi = diff(d$psqi), d_pvs = diff(d$total_pvs_volume),
             d_cso = diff(d$pvs_volume_centrumsemiovale),
             d_bg  = diff(d$pvs_volume_basalganglia))
}
block_permute <- function(x, block_size) {
  n <- length(x)
  x_ext <- c(x, x)
  starts <- sample(seq_len(n), size = ceiling(n / block_size), replace = TRUE)
  unlist(lapply(starts, function(s) x_ext[s:(s + block_size - 1)]))[seq_len(n)]
}
coupling <- function(x, y) {
  set.seed(42)
  r <- cor(x, y, method = "spearman")
  r_perm <- replicate(N_PERM, cor(block_permute(x, BLOCK_SIZE), y, method = "spearman"))
  list(r = r, p = (sum(abs(r_perm) >= abs(r)) + 1) / (N_PERM + 1))
}
cp_label <- function(cp) sprintf("italic(rho) == %.2f*','~~italic(p)[perm]~'%s'", cp$r, fmt_p(cp$p))

nat_delta <- make_delta(nat)

region_panel <- function(y_var, y_lab, cp) {
  d <- nat_delta
  fit <- lm(as.formula(paste(y_var, "~ d_psqi")), data = d)
  x_seq <- seq(min(d$d_psqi), max(d$d_psqi), length.out = 200)
  pr <- predict(fit, newdata = data.frame(d_psqi = x_seq), interval = "confidence")
  reg_df <- data.frame(x = x_seq, fit = pr[, "fit"], lo = pr[, "lwr"], hi = pr[, "upr"])
  y_rng <- range(c(d[[y_var]], reg_df$lo, reg_df$hi))
  y_top <- y_rng[2] + 0.2 * diff(y_rng)
  col <- "#4A8585"
  ggplot() +
    geom_hline(yintercept = 0, color = "grey40", linewidth = LW_THIN, linetype = "22") +
    geom_vline(xintercept = 0, color = "grey40", linewidth = LW_THIN, linetype = "22") +
    geom_ribbon(data = reg_df, aes(x = x, ymin = lo, ymax = hi), fill = col, alpha = 0.2) +
    geom_line(data = reg_df, aes(x = x, y = lo), color = col, linewidth = LW_THIN) +
    geom_line(data = reg_df, aes(x = x, y = hi), color = col, linewidth = LW_THIN) +
    geom_line(data = reg_df, aes(x = x, y = fit), color = col, linewidth = LW_MAIN) +
    geom_point(data = d, aes(x = d_psqi, y = .data[[y_var]], fill = phase, shape = phase),
               size = POINT_SIZE, stroke = POINT_STROKE, color = "black") +
    stat_layer(cp_label(cp), x = min(d$d_psqi), y = y_top) +
    scale_fill_manual(values = phase_colors, guide = "none") +
    scale_shape_manual(values = phase_shapes, guide = "none") +
    scale_x_continuous(breaks = pretty_breaks(n = 6)) +
    coord_cartesian(ylim = c(y_rng[1], y_top)) +
    labs(x = expression(Delta*"1-dPSQI"), y = y_lab) +
    theme_fig()
}

cp_cso <- coupling(nat_delta$d_psqi, nat_delta$d_cso)
cp_bg  <- coupling(nat_delta$d_psqi, nat_delta$d_bg)
s_region <- cowplot::plot_grid(
  label_panels(region_panel("d_cso", expression(Delta*"CSO-PVSV (mm"^3*")"), cp_cso),
               region_panel("d_bg",  expression(Delta*"BG-PVSV (mm"^3*")"),  cp_bg),
               ncol = 2, labels = c("a", "b")),
  phase_legend(phase_colors, phase_shapes, "Cycle phase"),
  ncol = 1, rel_heights = c(1, 0.1))
save_fig(s_region, "SuppFig5", 180, 80)

# ══════════════════════════════════════════════════════════════════════════════
# Spider plots: correlations with the 1-dPSQI components
# ══════════════════════════════════════════════════════════════════════════════
comp_vars <- c("C1_quality", "C2_latency", "C3_duration", "C4_efficiency",
               "C5_disturbances", "C7_daytime", "psqi")
comp_labs <- c("Subjective\nsleep quality", "Sleep\nlatency", "Sleep\nduration",
               "Habitual sleep\nefficiency", "Sleep\ndisturbances", "Daytime\ndysfunction",
               "Global\nscore")

# psqi_scored$session_num may use the session IDs of df or days within
# condition (1-30); both are mapped to the session IDs of df
map_components <- function(d, cond) {
  ps  <- psqi_scored[psqi_scored$condition == cond, ]
  ids <- ps$session_num
  if (!all(ids %in% d$session)) {
    if (!setequal(ids, seq_len(nrow(d)))) stop("psqi_scored: session_num cannot be matched for ", cond)
    ids <- d$session[order(d$session)][ids]
  }
  ps$session <- ids
  merge(d, ps[, c("session", setdiff(comp_vars, "psqi"))], by = "session")
}
nat_c <- map_components(nat, "natural")

spearman_row <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(r = unname(ct$estimate), p = ct$p.value)
}

# Radial scale from -r_max (centre) to +r_max (outer ring); r_max is shared
# within a figure so that panels are comparable
radial_max <- function(r_all) max(0.5, ceiling(max(abs(r_all)) * 2) / 2)

spider_panel <- function(r, p, color, tag, r_max = 0.5) {
  n   <- length(r)
  ang <- pi / 2 - seq(0, 2 * pi * (1 - 1 / n), length.out = n)   # clockwise from top
  rad <- function(v) (v + r_max) / (2 * r_max)                    # -r_max at centre
  step <- if (r_max <= 0.5) 0.25 else 0.5
  ring_vals <- seq(-r_max, r_max, by = step)
  lab_vals  <- ring_vals[seq(length(ring_vals), 2, by = -2)]       # every second ring
  ang_lab   <- pi / 2 - pi / n                                     # between first two axes
  rings <- do.call(rbind, lapply(ring_vals, function(v) {
    a <- seq(0, 2 * pi, length.out = 200)
    data.frame(x = rad(v) * cos(a), y = rad(v) * sin(a), ring = v)
  }))
  poly <- data.frame(x = rad(c(r, r[1])) * cos(c(ang, ang[1])),
                     y = rad(c(r, r[1])) * sin(c(ang, ang[1])))
  pts  <- data.frame(x = rad(r) * cos(ang), y = rad(r) * sin(ang), sig = p < 0.05)
  # Labels just outside the outer ring, aligned away from the plot centre
  labs_df <- data.frame(x = 1.1 * cos(ang), y = 1.1 * sin(ang),
                        lab = paste0(comp_labs, "\n", sprintf("%+.2f", r),
                                     ifelse(p < 0.05, "*", "")),
                        face = ifelse(p < 0.05, "bold", "plain"),
                        hjust = ifelse(cos(ang) > 0.2, 0, ifelse(cos(ang) < -0.2, 1, 0.5)),
                        vjust = ifelse(sin(ang) > 0.8, 0, ifelse(sin(ang) < -0.8, 1, 0.5)))

  ggplot() +
    geom_path(data = subset(rings, ring != 0), aes(x = x, y = y, group = ring),
              color = "grey70", linewidth = LW_THIN) +
    geom_path(data = subset(rings, ring == 0), aes(x = x, y = y),
              color = "grey40", linewidth = LW_THIN, linetype = "22") +
    geom_segment(data = data.frame(xend = cos(ang), yend = sin(ang)),
                 aes(x = 0, y = 0, xend = xend, yend = yend),
                 color = "grey70", linewidth = LW_THIN) +
    geom_polygon(data = poly, aes(x = x, y = y), fill = color, alpha = 0.2,
                 color = color, linewidth = LW_MAIN) +
    geom_point(data = pts, aes(x = x, y = y, fill = sig), shape = 21,
               size = 2, stroke = 0.6, color = color) +
    geom_text(data = labs_df, aes(x = x, y = y, label = lab, fontface = face,
                                  hjust = hjust, vjust = vjust),
              size = pt_to_size(TXT), lineheight = 0.9, family = FONT_FAMILY) +
    annotate("text", x = rad(lab_vals) * cos(ang_lab), y = rad(lab_vals) * sin(ang_lab),
             label = ifelse(lab_vals == 0, "0", sprintf("%+.1f", lab_vals)),
             hjust = -0.15, vjust = -0.2,
             size = pt_to_size(TXT - 1), family = FONT_FAMILY, color = "grey25") +
    annotate("text", x = -2.1, y = 1.95, label = tag, hjust = 0, vjust = 1,
             size = pt_to_size(TXT), family = FONT_FAMILY) +
    scale_fill_manual(values = c(`TRUE` = unname(color), `FALSE` = "white"), guide = "none") +
    coord_equal(xlim = c(-2.1, 2.1), ylim = c(-1.75, 1.95), clip = "off") +
    theme_void(base_family = FONT_FAMILY) +
    theme(plot.margin = margin(12, 4, 4, 4),
          plot.background = element_rect(fill = "white", color = NA))
}

spider_legend <- extract_legend(
  ggplot(data.frame(x = 1:2, y = 1, s = factor(c("p < 0.05", "n.s."), levels = c("p < 0.05", "n.s."))),
         aes(x, y, fill = s)) +
    geom_point(shape = 21, size = 2, stroke = 0.6, color = "black") +
    scale_fill_manual(values = c("p < 0.05" = "grey30", "n.s." = "white"),
                      name = "Spearman correlation") +
    theme_fig() + theme(legend.position = "bottom"))

# Hormones x components (physiological cycle)
hormone_spiders <- function(d, name) {
  sts <- lapply(hormones, function(h) t(sapply(comp_vars, function(v) spearman_row(d[[h]], d[[v]]))))
  r_max <- radial_max(unlist(lapply(sts, function(st) st[, "r"])))
  panels <- lapply(seq_along(hormones), function(k)
    spider_panel(sts[[k]][, "r"], sts[[k]][, "p"], hormone_colors[horm_lab[hormones[k]]],
                 horm_lab[hormones[k]], r_max))
  fig <- cowplot::plot_grid(do.call(label_panels, c(panels, list(ncol = 2, labels = letters[1:4]))),
                            spider_legend, ncol = 1, rel_heights = c(1, 0.05))
  save_fig(fig, name, 180, 185)
}
hormone_spiders(nat_c, "SuppFig3")

# PVSV x components, both conditions
pvs_component_stats <- function(d, outlier = numeric(0)) {
  if (PVS_SPIDER_DELTA) {
    d <- d[order(d$session), ]
    dd <- data.frame(session = d$session[-1], d_pvs = diff(d$total_pvs_volume))
    for (v in comp_vars) dd[[v]] <- diff(d[[v]])
    dd <- dd[!(dd$session %in% c(outlier, outlier + 1)), ]
    t(sapply(comp_vars, function(v) spearman_row(dd$d_pvs, dd[[v]])))
  } else {
    d <- d[!(d$session %in% outlier), ]
    t(sapply(comp_vars, function(v) spearman_row(d$total_pvs_volume, d[[v]])))
  }
}
oc_c_all <- map_components(oc, "OC")
st_nat <- pvs_component_stats(nat_c)
st_oc  <- pvs_component_stats(oc_c_all, outlier_oc)
s_pvs <- cowplot::plot_grid(
  label_panels(spider_panel(st_nat[, "r"], st_nat[, "p"], "#4A8585", "Physiological menstrual cycle",
                            radial_max(c(st_nat[, "r"], st_oc[, "r"]))),
               spider_panel(st_oc[, "r"],  st_oc[, "p"],  "#4A8585", "OCP",
                            radial_max(c(st_nat[, "r"], st_oc[, "r"]))),
               ncol = 2, labels = c("a", "b")),
  spider_legend, ncol = 1, rel_heights = c(1, 0.08))
save_fig(s_pvs, "SuppFig4", 180, 100)
