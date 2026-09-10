# Figure: trajectory + conditional slopes + interaction plots (5 panels)

# ---- DATA PREPARATION ----
nat <- df[df$condition == "natural", ]
nat <- nat[order(nat$session), ]

nat$cycle_day <- ifelse(nat$session >= 11,
                        nat$session - 10,
                        nat$session + 20)

nat$phase <- cut(nat$cycle_day,
                 breaks = c(0, 4, 11, 14, 30),
                 labels = c("Menses","Follicular",
                            "Periovulatory","Luteal"))

nat <- nat[order(nat$cycle_day), ]

nat$log_E2  <- log(nat$`Estradiol (LCMS; pg/mL)`    + 1)
nat$log_P4  <- log(nat$`Progesterone (LCMS; ng/mL)`  + 1)
nat$log_LH  <- log(nat$`LH (mIU/ml)`  + 1)
nat$log_FSH <- log(nat$`FSH (mIU/ml)` + 1)

nat_ar <- nat %>%
  arrange(cycle_day) %>%
  mutate(time = seq_len(n()))

# ---- INTERACTION MODELS ----
library(nlme)
library(lmtest)

dw_lh <- dwtest(lm(total_pvs_volume ~ psqi * log_LH, data = nat_ar),
                alternative = "two.sided")
if (dw_lh$p.value < 0.05) {
  m3_ar_lh <- gls(total_pvs_volume ~ psqi * log_LH,
                  data        = nat_ar,
                  correlation = corAR1(form = ~time),
                  method      = "ML")
} else {
  m3_ar_lh <- lm(total_pvs_volume ~ psqi * log_LH, data = nat_ar)
}

dw_fsh <- dwtest(lm(total_pvs_volume ~ psqi * log_FSH, data = nat_ar),
                 alternative = "two.sided")
if (dw_fsh$p.value < 0.05) {
  m4_ar_fsh <- gls(total_pvs_volume ~ psqi * log_FSH,
                   data        = nat_ar,
                   correlation = corAR1(form = ~time),
                   method      = "ML")
} else {
  m4_ar_fsh <- lm(total_pvs_volume ~ psqi * log_FSH, data = nat_ar)
}

library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
library(ggthemes)

# Unicode rendering for PDF export
if (!requireNamespace("showtext", quietly = TRUE)) {
  install.packages("showtext")
}
library(showtext)

save_pdf <- function(filename, plot, width, height, units = "in") {
  showtext_opts(dpi = 72)
  showtext_auto()
  ggsave(filename, plot, width = width, height = height,
         units = units, bg = "white")
  showtext_auto(FALSE)
}

# For italic/subscript stat annotations
if (!requireNamespace("ggtext", quietly = TRUE)) {
  install.packages("ggtext")
}
library(ggtext)

# ---- COLOURS ----
hormone_colors <- c("LH" = "#BF8A3D", "FSH" = "#4F9483")
# Colour-blind-safe palette
phase_colors_nat <- c(
  "Menses"        = "#D55E00",
  "Follicular"    = "#0072B2",
  "Periovulatory" = "#009E73",
  "Luteal"        = "#CC79A7"
)
reg_colors <- c("lh" = "#BF8A3D", "fsh" = "#4F9483")
trajectory_line_colors <- c("LH" = "#A6762C", "FSH" = "#3D7A67")

# Base font size (must match delta_reg_plot.R)
TITLE_BASE_SIZE <- 16
LEGEND_TEXT_BASE_SIZE <- 13 * (TITLE_BASE_SIZE / 14)

theme_nature <- function(base_size = TITLE_BASE_SIZE) {
  theme_foundation(base_size = base_size, base_family = "sans") +
    theme(
      text              = element_text(family = "sans"),
      plot.title        = element_text(face = "bold", size = base_size + 4,
                                       hjust = 0, margin = margin(b = 2),
                                       color = "grey10", lineheight = 1.2),
      axis.title        = element_text(size = base_size, color = "black"),
      axis.text         = element_text(size = base_size - 2, color = "black"),
      axis.ticks        = element_line(color = "black", linewidth = 0.25),
      axis.ticks.length = unit(2.5, "pt"),
      axis.line         = element_line(color = "black", linewidth = 0.35),
      panel.background  = element_rect(fill = "white", color = NA),
      panel.border      = element_blank(),
      panel.grid        = element_blank(),
      legend.position   = "right",
      legend.key.size   = unit(12, "pt"),
      legend.key.width  = unit(28, "pt"),
      legend.text       = element_text(size = LEGEND_TEXT_BASE_SIZE),
      legend.title      = element_text(size = TITLE_BASE_SIZE, face = "bold"),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.spacing.y  = unit(1, "pt"),
      plot.background   = element_rect(fill = "white", color = NA),
      plot.margin       = margin(8, 14, 6, 6)
    )
}

# Text scaling, referenced to delta_reg_plot.R Figure 1
ref_panel_w <- 18   / (1 + 1 + 0.22)
ref_panel_h <- 13.0 / (1 + 1 + 0.12)
ref_panel_area <- ref_panel_w * ref_panel_h

panel_w_wide <- 14   / (1 + 0.14)
panel_h_wide <- 15.9 * (0.36 / 1.87)

panel_w_int <- 14   / (1 + 1 + 0.28)
panel_h_int <- 15.9 * (0.72 / 1.87)

# Corrects for different ggsave() canvas widths between scripts
REF_CANVAS_WIDTH_IN <- 18
OWN_CANVAS_WIDTH_IN <- 14
CANVAS_CORR <- OWN_CANVAS_WIDTH_IN / REF_CANVAS_WIDTH_IN

text_scale_wide <- sqrt((panel_w_wide * panel_h_wide) / ref_panel_area) * CANVAS_CORR
text_scale_int  <- sqrt((panel_w_int  * panel_h_int ) / ref_panel_area) * CANVAS_CORR

# Minimum readable text size
TEXT_SCALE_FLOOR <- 0.87
text_scale_wide <- pmax(text_scale_wide, TEXT_SCALE_FLOOR)
text_scale_int  <- pmax(text_scale_int,  TEXT_SCALE_FLOOR)

STAT_SIZE_FRACTION_OF_TITLE <- 0.90
BASE_STAT_SIZE  <- (TITLE_BASE_SIZE + 4) * STAT_SIZE_FRACTION_OF_TITLE / 2.845276
BASE_LABEL_SIZE <- 5.2
BASE_POINT_SIZE   <- 3.5
BASE_POINT_STROKE <- 0.5

# Y-axis scale matches delta_reg_plot.R
panel_h_fig2_main <- 7.0 / (1 + 0.08)

ref_pvs_lim  <- c(floor(min(nat$total_pvs_volume)/50)*50,
                  ceiling(max(nat$total_pvs_volume)/50)*50)
ref_y_expand <- 0.18
ref_pvs_span <- diff(ref_pvs_lim) * (1 + ref_y_expand)
PVSV_UNITS_PER_INCH <- ref_pvs_span / panel_h_fig2_main

mod_pvs_span   <- PVSV_UNITS_PER_INCH * panel_h_int
mod_pvs_center <- mean(range(nat$total_pvs_volume, na.rm = TRUE))
pvs_lim <- c(mod_pvs_center - mod_pvs_span/2, mod_pvs_center + mod_pvs_span/2)

phases_nat <- data.frame(
  xmin  = c(1,  12),
  xmax  = c(4,  14),
  label = c("M", "O")
)

nat_ts <- nat_ar %>%
  arrange(cycle_day) %>%
  mutate(
    lh_z  = as.numeric(scale(log_LH)),
    fsh_z = as.numeric(scale(log_FSH))
  )

smooth_series <- function(x, y, n = 300) {
  keep <- !is.na(x) & !is.na(y)
  s <- spline(x = x[keep], y = y[keep], n = n)
  data.frame(x = s$x, y = s$y)
}

# Conditional slope coefficients
b_psqi_lh  <- coef(m3_ar_lh)["psqi"]
b_inter_lh <- coef(m3_ar_lh)["psqi:log_LH"]
vc_lh      <- vcov(m3_ar_lh)

b_psqi_fsh  <- coef(m4_ar_fsh)["psqi"]
b_inter_fsh <- coef(m4_ar_fsh)["psqi:log_FSH"]
vc_fsh      <- vcov(m4_ar_fsh)

slope_df <- nat_ar %>%
  arrange(cycle_day) %>%
  mutate(
    cond_slope_lh = b_psqi_lh + b_inter_lh * log_LH,
    se_lh = sqrt(pmax(0,
                      vc_lh["psqi","psqi"] +
                        log_LH^2 * vc_lh["psqi:log_LH","psqi:log_LH"] +
                        2 * log_LH * vc_lh["psqi","psqi:log_LH"])),
    ci_lo_lh = cond_slope_lh - 1.96 * se_lh,
    ci_hi_lh = cond_slope_lh + 1.96 * se_lh,
    sig_lh   = !is.na(ci_lo_lh) & !is.na(ci_hi_lh) &
      (ci_lo_lh > 0 | ci_hi_lh < 0),
    
    cond_slope_fsh = b_psqi_fsh + b_inter_fsh * log_FSH,
    se_fsh = sqrt(pmax(0,
                       vc_fsh["psqi","psqi"] +
                         log_FSH^2 * vc_fsh["psqi:log_FSH","psqi:log_FSH"] +
                         2 * log_FSH * vc_fsh["psqi","psqi:log_FSH"])),
    ci_lo_fsh = cond_slope_fsh - 1.96 * se_fsh,
    ci_hi_fsh = cond_slope_fsh + 1.96 * se_fsh,
    sig_fsh   = !is.na(ci_lo_fsh) & !is.na(ci_hi_fsh) &
      (ci_lo_fsh > 0 | ci_hi_fsh < 0)
  )

shared_ylim <- range(c(slope_df$ci_lo_lh, slope_df$ci_hi_lh,
                       slope_df$ci_lo_fsh, slope_df$ci_hi_fsh), na.rm = TRUE)
shared_ylim <- shared_ylim + diff(shared_ylim) * c(-0.08, 0.05)

# ---- PANEL A: TRAJECTORY ----
smooth_fsh_a <- smooth_series(nat_ts$cycle_day, nat_ts$fsh_z, n = 300)
smooth_lh_a  <- smooth_series(nat_ts$cycle_day, nat_ts$lh_z,  n = 300)

panel_a <- ggplot() +
  geom_rect(data = phases_nat,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey93", alpha = 1, inherit.aes = FALSE) +
  geom_text(data = phases_nat,
            aes(x = (xmin + xmax) / 2, y = 4.55, label = label),
            size = BASE_LABEL_SIZE * text_scale_wide,
            color = "black", fontface = "bold", inherit.aes = FALSE) +
  geom_ribbon(data = smooth_fsh_a, aes(x = x, ymin = 0, ymax = y),
              fill = hormone_colors["FSH"], color = NA, alpha = 0.18,
              inherit.aes = FALSE) +
  geom_ribbon(data = smooth_lh_a, aes(x = x, ymin = 0, ymax = y),
              fill = hormone_colors["LH"], color = NA, alpha = 0.20,
              inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey65", linewidth = 0.22) +
  geom_line(data = smooth_fsh_a, aes(x = x, y = y),
            color = trajectory_line_colors["FSH"], linewidth = 1.0,
            inherit.aes = FALSE) +
  geom_line(data = smooth_lh_a, aes(x = x, y = y),
            color = trajectory_line_colors["LH"], linewidth = 1.0,
            inherit.aes = FALSE) +
  scale_x_continuous(breaks = 1:30,
                     limits = c(1, 30), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-2, 4, by = 1)) +
  coord_cartesian(ylim = c(-2.0, 4.8), clip = "off") +
  labs(title = "A   FSH and LH across the physiological menstrual cycle",
       x = "Cycle day", y = "z-score (log-transformed)") +
  theme_nature(base_size = TITLE_BASE_SIZE * text_scale_wide) +
  theme(plot.margin  = margin(8, 6, 2, 6))

leg_a <- ggdraw() +
  draw_text("Gonadotropins",
            x = 0.05, y = 0.78, hjust = 0, vjust = 1,
            size = LEGEND_TEXT_BASE_SIZE * text_scale_wide, fontface = "bold",
            color = "grey15", lineheight = 0.9) +
  draw_grob(grid::rectGrob(
    x = 0.08, y = 0.58, width = 0.22, height = 0.10,
    just = c("left","centre"),
    gp = grid::gpar(fill = scales::alpha(hormone_colors["FSH"], 0.45), col = NA)
  )) +
  draw_text("FSH", x = 0.35, y = 0.58, hjust = 0, vjust = 0.5,
            size = LEGEND_TEXT_BASE_SIZE * text_scale_wide, color = "grey15") +
  draw_grob(grid::rectGrob(
    x = 0.08, y = 0.44, width = 0.22, height = 0.10,
    just = c("left","centre"),
    gp = grid::gpar(fill = scales::alpha(hormone_colors["LH"], 0.55), col = NA)
  )) +
  draw_text("LH",  x = 0.35, y = 0.44, hjust = 0, vjust = 0.5,
            size = LEGEND_TEXT_BASE_SIZE * text_scale_wide, color = "grey15")

# ---- PANELS D/E: INTERACTION (FSH, LH) ----
make_interaction <- function(data, x_var, y_var, mod_var,
                             x_label, y_label, mod_label,
                             title, mod_color, phase_var, phase_colors,
                             y_limits = NULL, x_limits = NULL, text_scale = 1) {
  
  keep       <- complete.cases(data[[x_var]], data[[y_var]], data[[mod_var]])
  data_clean <- data[keep, ]
  
  mod_mean <- mean(data_clean[[mod_var]], na.rm = TRUE)
  mod_sd   <- sd(data_clean[[mod_var]],   na.rm = TRUE)
  
  x_seq <- seq(min(data_clean[[x_var]], na.rm = TRUE),
               max(data_clean[[x_var]], na.rm = TRUE), length.out = 200)
  
  fit_model <- lm(
    as.formula(paste0("`", y_var, "` ~ `", x_var, "` * `", mod_var, "`")),
    data = data_clean
  )
  
  mod_levels <- list(
    list(label = "+1 SD", val = mod_mean + mod_sd),
    list(label = "Mean",  val = mod_mean),
    list(label = "-1 SD", val = mod_mean - mod_sd)
  )
  
  pred_lines <- lapply(mod_levels, function(m) {
    nd   <- setNames(data.frame(x_seq, m$val), c(x_var, mod_var))
    pred <- predict(fit_model, newdata = nd, interval = "confidence", level = 0.95)
    df   <- data.frame(x = x_seq, fit = pred[,"fit"],
                       lower = pred[,"lwr"], upper = pred[,"upr"],
                       group = m$label)
    df[complete.cases(df), ]
  })
  
  pred_df       <- do.call(rbind, pred_lines)
  pred_df$group <- factor(pred_df$group, levels = c("+1 SD","Mean","-1 SD"))
  pred_df       <- pred_df %>% arrange(group, x)
  
  line_types  <- c("+1 SD" = "solid",  "Mean" = "dashed",  "-1 SD" = "dotted")
  line_widths <- c("+1 SD" = 1.0,      "Mean" = 0.8,       "-1 SD" = 0.8)
  
  band_df <- pred_df %>%
    group_by(x) %>%
    summarise(lower = min(lower), upper = max(upper), .groups = "drop") %>%
    arrange(x)
  
  # Cap the ribbon/fit lines below the stat box instead of extending the axis
  if (!is.null(y_limits)) {
    ribbon_ceiling  <- y_limits[2] - diff(y_limits) * 0.08
    band_df$upper   <- pmin(band_df$upper, ribbon_ceiling)
    pred_df$fit     <- pmin(pred_df$fit,   ribbon_ceiling)
  }
  
  p <- ggplot() +
    geom_ribbon(data = band_df,
                aes(x = x, ymin = lower, ymax = upper),
                fill = mod_color, alpha = 0.15, color = NA) +
    geom_line(data = pred_df,
              aes(x = x, y = fit, linetype = group, linewidth = group),
              color = mod_color) +
    geom_point(data = data_clean,
               aes(x = .data[[x_var]], y = .data[[y_var]],
                   fill = .data[[phase_var]]),
               size = BASE_POINT_SIZE * text_scale, alpha = 0.88, shape = 21,
               stroke = BASE_POINT_STROKE * text_scale, color = "white") +
    scale_fill_manual(values = phase_colors, name = NULL, guide = "none") +
    scale_linetype_manual(
      values = line_types,
      name   = paste0(mod_label, " (log, \u00b11 SD)"),
      guide  = guide_legend(override.aes = list(linewidth = c(1.0, 0.8, 0.8)))
    ) +
    scale_linewidth_manual(
      values = line_widths,
      name   = paste0(mod_label, " (log, \u00b11 SD)")
    ) +
    labs(title = title, x = x_label, y = y_label) +
    theme_nature(base_size = TITLE_BASE_SIZE * text_scale) +
    theme(legend.key.width = unit(28, "pt"))
  
  if (!is.null(x_limits)) {
    p <- p + scale_x_continuous(limits = x_limits, oob = squish,
                                breaks = pretty_breaks(n = 4))
  } else {
    p <- p + scale_x_continuous(breaks = pretty_breaks(n = 4))
  }
  
  if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(limits = y_limits, oob = squish,
                                breaks = pretty_breaks(n = 4))
  } else {
    p <- p + scale_y_continuous(breaks = pretty_breaks(n = 4))
  }
  return(p)
}

get_inter_stats <- function(model, inter_term) {
  if (inherits(model, "gls")) {
    s   <- summary(model)$tTable
    b   <- s[inter_term, "Value"]
    se  <- s[inter_term, "Std.Error"]
    p   <- s[inter_term, "p-value"]
    df_ <- s[inter_term, "DF"]
    ci  <- b + c(-1,1) * qt(0.975, df_) * se
  } else {
    s  <- summary(model)$coefficients
    b  <- s[inter_term, "Estimate"]
    se <- s[inter_term, "Std. Error"]
    p  <- s[inter_term, "Pr(>|t|)"]
    ci <- confint(model)[inter_term, ]
  }
  list(beta=b, se=se, p=p, ci_lo=ci[1], ci_hi=ci[2])
}

p_inter_all <- sapply(
  list(lh=m3_ar_lh, fsh=m4_ar_fsh,
       e2=lm(total_pvs_volume ~ psqi * log_E2, data=nat_ar),
       p4=lm(total_pvs_volume ~ psqi * log_P4, data=nat_ar)),
  function(m) {
    nm <- grep("psqi:", names(coef(m)), value=TRUE)
    if (inherits(m,"gls")) summary(m)$tTable[nm,"p-value"]
    else summary(m)$coefficients[nm,"Pr(>|t|)"]
  }
)
p_inter_fdr <- p.adjust(p_inter_all, method="BH")

cf_lh  <- get_inter_stats(m3_ar_lh,  "psqi:log_LH")
cf_fsh <- get_inter_stats(m4_ar_fsh, "psqi:log_FSH")

# Number formatting: normal decimal point, 2 sig. digits for p-values
mid_dot <- function(x) x
fmt2 <- function(x) mid_dot(sprintf("%.2f", x))
fmt_p2sig <- function(p) {
  if (p < 0.0001) return("<0.0001")
  rounded <- signif(p, 2)
  if (rounded <= 0) return("<0.0001")
  decimals <- max(0, -floor(log10(rounded)) + 1)
  mid_dot(sprintf(paste0("%.", decimals, "f"), rounded))
}

sig_inter <- function(p) ifelse(p<0.001,"\\*\\*\\*",ifelse(p<0.01,"\\*\\*",
                                                           ifelse(p<0.05,"\\*",ifelse(p<0.10,"\u2020",""))))

mid_dot_raw   <- function(x) x
fmt2_raw      <- function(x) mid_dot_raw(sprintf("%.2f", x))
sig_inter_raw <- function(p) ifelse(p<0.001,"***",ifelse(p<0.01,"**",
                                                         ifelse(p<0.05,"*",ifelse(p<0.10,"\u2020",""))))
fmt_p2sig_raw <- function(p) {
  if (p < 0.0001) return("<0.0001")
  rounded <- signif(p, 2)
  if (rounded <= 0) return("<0.0001")
  decimals <- max(0, -floor(log10(rounded)) + 1)
  mid_dot_raw(sprintf(paste0("%.", decimals, "f"), rounded))
}

# Multi-line stat box (grid::textGrob-based, see delta_reg_plot.R)
build_rich_label_grob_ml <- function(lines, x=unit(0.97,"npc"), y=unit(0.93,"npc"),
                                     hjust=c("right","left"),
                                     size=14, sub_scale=0.70, sub_dy_pt=3,
                                     line_gap_pt=NULL, fontface="bold", color="grey20",
                                     fill=scales::alpha("white", 0.75), pad_pt=4) {
  hjust <- match.arg(hjust)
  if (is.null(line_gap_pt)) line_gap_pt <- size * 1.35
  
  line_grobs  <- vector("list", length(lines))
  line_widths <- vector("list", length(lines))
  for (li in seq_along(lines)) {
    segs <- lines[[li]]
    gps <- lapply(segs, function(seg) {
      fs <- if (isTRUE(seg$sub)) size * sub_scale else size
      ff <- if (isTRUE(seg$italic)) "bold.italic" else fontface
      grid::gpar(fontsize=fs, fontface=ff, col=color)
    })
    measure_grobs <- Map(function(seg, gp) grid::textGrob(seg$text, gp=gp), segs, gps)
    widths        <- lapply(measure_grobs, function(g) grid::unit(1, "grobwidth", g))
    total_width   <- Reduce(`+`, widths)
    line_widths[[li]] <- total_width
    
    this_y <- y - grid::unit((li-1) * line_gap_pt, "pt")
    start_x <- if (hjust == "right") x - total_width else x
    
    cur_x <- start_x
    segs_grobs <- vector("list", length(segs))
    for (i in seq_along(segs)) {
      seg <- segs[[i]]
      seg_y <- if (isTRUE(seg$sub)) this_y - grid::unit(sub_dy_pt, "pt") else this_y
      segs_grobs[[i]] <- grid::textGrob(segs[[i]]$text, x=cur_x, y=seg_y,
                                        hjust=0, vjust=1, gp=gps[[i]])
      cur_x <- cur_x + widths[[i]]
    }
    line_grobs[[li]] <- do.call(grid::grobTree, segs_grobs)
  }
  
  do.call(grid::grobTree, line_grobs)
}

lines_inter_lh <- list(
  list(list(text="b", sub=FALSE, italic=TRUE),
       list(text=paste0("=", fmt2_raw(cf_lh$beta),
                        ", 95% (CI ", fmt2_raw(cf_lh$ci_lo), ", ",
                        fmt2_raw(cf_lh$ci_hi), "), SE=", fmt2_raw(cf_lh$se), ","), sub=FALSE)),
  list(list(text="p", sub=FALSE, italic=TRUE),
       list(text=paste0("=", fmt_p2sig_raw(cf_lh$p), sig_inter_raw(cf_lh$p),
                        ", "), sub=FALSE),
       list(text="p", sub=FALSE, italic=TRUE),
       list(text="FDR", sub=TRUE),
       list(text=paste0("=", fmt_p2sig_raw(p_inter_fdr["lh"]), sig_inter_raw(p_inter_fdr["lh"])), sub=FALSE))
)
lines_inter_fsh <- list(
  list(list(text="b", sub=FALSE, italic=TRUE),
       list(text=paste0("=", fmt2_raw(cf_fsh$beta),
                        ", 95% (CI ", fmt2_raw(cf_fsh$ci_lo), ", ",
                        fmt2_raw(cf_fsh$ci_hi), "), SE=", fmt2_raw(cf_fsh$se), ","), sub=FALSE)),
  list(list(text="p", sub=FALSE, italic=TRUE),
       list(text=paste0("=", fmt_p2sig_raw(cf_fsh$p), sig_inter_raw(cf_fsh$p),
                        ", "), sub=FALSE),
       list(text="p", sub=FALSE, italic=TRUE),
       list(text="FDR", sub=TRUE),
       list(text=paste0("=", fmt_p2sig_raw(p_inter_fdr["fsh"]), sig_inter_raw(p_inter_fdr["fsh"])), sub=FALSE))
)

panel_b_base <- make_interaction(
  nat_ar, "psqi", "total_pvs_volume", "log_FSH",
  "1-dPSQI", "PVSV (mm\u00b3)", "FSH",
  "D   PVSV ~ 1-dPSQI \u00d7 FSH",
  reg_colors["fsh"], "phase", phase_colors_nat,
  y_limits = pvs_lim,
  text_scale = text_scale_int
) + annotation_custom(
  build_rich_label_grob_ml(lines_inter_fsh, y=unit(0.97,"npc"),
                           size=BASE_STAT_SIZE * text_scale_int * 2.845276),
  xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

panel_c_base <- make_interaction(
  nat_ar, "psqi", "total_pvs_volume", "log_LH",
  "1-dPSQI", "PVSV (mm\u00b3)", "LH",
  "E   PVSV ~ 1-dPSQI \u00d7 LH",
  reg_colors["lh"], "phase", phase_colors_nat,
  y_limits = pvs_lim,
  text_scale = text_scale_int
) + annotation_custom(
  build_rich_label_grob_ml(lines_inter_lh, y=unit(0.97,"npc"),
                           size=BASE_STAT_SIZE * text_scale_int * 2.845276),
  xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

panel_b <- panel_b_base + theme(legend.position = "none")
panel_c <- panel_c_base + theme(legend.position = "none")

leg_d <- cowplot::get_legend(
  panel_b_base +
    theme(legend.position  = "right",
          legend.key.width = unit(36, "pt"),
          legend.text      = element_text(size = LEGEND_TEXT_BASE_SIZE * text_scale_int),
          legend.title     = element_text(size = TITLE_BASE_SIZE * text_scale_int, face = "bold"),
          legend.box.margin = margin(0, 0, 0, 8))
)
leg_e <- cowplot::get_legend(
  panel_c_base +
    theme(legend.position  = "right",
          legend.key.width = unit(36, "pt"),
          legend.text      = element_text(size = LEGEND_TEXT_BASE_SIZE * text_scale_int),
          legend.title     = element_text(size = TITLE_BASE_SIZE * text_scale_int, face = "bold"),
          legend.box.margin = margin(0, 0, 0, 8))
)
leg_de <- cowplot::plot_grid(leg_d, leg_e, ncol = 1,
                             align = "v", rel_heights = c(1, 1))

# ---- PANELS B/C: CONDITIONAL SLOPE (FSH, LH) ----
make_slope_panel <- function(data, slope_var, ci_lo_var, ci_hi_var,
                             hormone_color, title_label, spline_n = 300,
                             text_scale = 1) {
  
  data_ord <- data[order(data$cycle_day), ]
  data_ord <- data_ord[!is.na(data_ord[[slope_var]]) &
                         !is.na(data_ord[[ci_lo_var]]) &
                         !is.na(data_ord[[ci_hi_var]]), ]
  
  spline_slope <- spline(x = data_ord$cycle_day, y = data_ord[[slope_var]],
                         n = spline_n)
  spline_lo    <- spline(x = data_ord$cycle_day, y = data_ord[[ci_lo_var]],
                         n = spline_n)
  spline_hi    <- spline(x = data_ord$cycle_day, y = data_ord[[ci_hi_var]],
                         n = spline_n)
  
  smooth_df <- data.frame(
    cycle_day = spline_slope$x,
    slope     = spline_slope$y,
    ci_lo     = pmax(spline_lo$y, shared_ylim[1]),
    ci_hi     = spline_hi$y
  )
  
  ggplot() +
    annotate("rect", xmin = 1, xmax = 4,
             ymin = -Inf, ymax = Inf, fill = "grey93", alpha = 1) +
    annotate("text", x = 2.5, y = Inf, label = "M",
             vjust = 1.5, size = BASE_LABEL_SIZE * text_scale,
             fontface = "bold", color = "black") +
    annotate("rect", xmin = 12, xmax = 14,
             ymin = -Inf, ymax = Inf, fill = "grey93", alpha = 1) +
    annotate("text", x = 13, y = Inf, label = "O",
             vjust = 1.5, size = BASE_LABEL_SIZE * text_scale,
             fontface = "bold", color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey55", linewidth = 0.30) +
    geom_ribbon(data = smooth_df,
                aes(x = cycle_day, ymin = ci_lo, ymax = ci_hi),
                fill = hormone_color, alpha = 0.20, color = NA,
                inherit.aes = FALSE) +
    geom_line(data = smooth_df,
              aes(x = cycle_day, y = slope),
              color = hormone_color, linewidth = 1.0, inherit.aes = FALSE) +
    geom_point(data = data,
               aes(x = cycle_day, y = .data[[slope_var]], fill = phase),
               size = BASE_POINT_SIZE * text_scale, shape = 21,
               stroke = BASE_POINT_STROKE * text_scale,
               color = "white", alpha = 0.92, inherit.aes = FALSE) +
    scale_fill_manual(values = phase_colors_nat, name = "Cycle phase") +
    scale_x_continuous(breaks = 1:30,
                       limits = c(1, 30), expand = c(0, 0)) +
    scale_y_continuous(breaks = pretty_breaks(n = 3)) +
    coord_cartesian(ylim = shared_ylim, clip = "on") +
    labs(title = title_label, x = "Cycle day",
         y = "Slope") +
    theme_nature(base_size = TITLE_BASE_SIZE * text_scale) +
    theme(legend.position = "none",
          plot.margin     = margin(14, 6, 6, 6))
}

panel_d <- make_slope_panel(
  slope_df, "cond_slope_fsh", "ci_lo_fsh", "ci_hi_fsh",
  hormone_colors["FSH"],
  "B   Simple slope: 1-dPSQI \u2192 PVSV  [FSH \u00b7 OLS]",
  text_scale = text_scale_wide
)
panel_e <- make_slope_panel(
  slope_df, "cond_slope_lh", "ci_lo_lh", "ci_hi_lh",
  hormone_colors["LH"],
  "C   Simple slope: 1-dPSQI \u2192 PVSV  [LH \u00b7 OLS]",
  text_scale = text_scale_wide
)

make_phase_legend_horiz <- function(phase_colors, text_scale = 1) {
  df <- data.frame(
    phase = factor(names(phase_colors), levels = names(phase_colors)),
    x = seq_along(phase_colors), y = 1
  )
  p <- ggplot(df, aes(x = x, y = y, color = phase)) +
    geom_point(size = 5.0 * text_scale, shape = 16) +
    scale_color_manual(values = phase_colors, name = NULL) +
    guides(color = guide_legend(
      nrow = 1,
      override.aes = list(size = 5.0 * text_scale),
      label.position = "right",
      keywidth  = unit(10, "pt"),
      keyheight = unit(10, "pt")
    )) +
    theme_void() +
    theme(
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.text       = element_text(size = 16 * (TITLE_BASE_SIZE/14) * text_scale, color = "grey15"),
      legend.key        = element_blank(),
      legend.spacing.x  = unit(3, "pt"),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.background   = element_rect(fill = "white", color = NA)
    )
  cowplot::get_legend(p)
}

leg_phase <- make_phase_legend_horiz(phase_colors_nat, text_scale = text_scale_wide)

# ---- ASSEMBLE ----
aligned <- cowplot::align_plots(
  panel_a, panel_b, panel_c, panel_d, panel_e,
  align = "v", axis = "l"
)

panel_a_al <- aligned[[1]]
panel_b_al <- aligned[[2]]
panel_c_al <- aligned[[3]]
panel_d_al <- aligned[[4]]
panel_e_al <- aligned[[5]]

# Second pass: align A/D/E panel heights
aligned_abc <- cowplot::align_plots(
  panel_a_al, panel_d_al, panel_e_al,
  align = "h", axis = "tb"
)
panel_a_al <- aligned_abc[[1]]
panel_d_al <- aligned_abc[[2]]
panel_e_al <- aligned_abc[[3]]

row1 <- cowplot::plot_grid(
  panel_a_al, leg_a,
  ncol = 2, rel_widths = c(1, 0.14), align = "h", axis = "tb"
)

spacer_align <- ggplot() + theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

row2 <- cowplot::plot_grid(
  panel_d_al, spacer_align,
  ncol = 2, rel_widths = c(1, 0.14), align = "h", axis = "tb"
)
row3 <- cowplot::plot_grid(
  panel_e_al, spacer_align,
  ncol = 2, rel_widths = c(1, 0.14), align = "h", axis = "tb"
)
row4 <- cowplot::plot_grid(
  panel_b_al, panel_c_al, leg_de,
  ncol = 3, rel_widths = c(1, 1, 0.28), align = "h", axis = "tb"
)

figure1 <- cowplot::plot_grid(
  row1, row2, row3, row4, leg_phase,
  ncol        = 1,
  rel_heights = c(0.36, 0.36, 0.36, 0.72, 0.07)
)

ggsave("figure1.png", figure1, width = 14, height = 15.9, dpi = 600, bg = "white")
save_pdf("figure1.pdf", figure1, width = 14, height = 15.9)