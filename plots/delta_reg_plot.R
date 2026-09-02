library(ggplot2)
library(dplyr)
library(ggthemes)
library(gridExtra)
library(grid)
library(nlme)
library(lmtest)

# Unicode rendering for PDF export
if (!requireNamespace("showtext", quietly = TRUE)) install.packages("showtext")
library(showtext)
if (!requireNamespace("ggtext", quietly = TRUE)) install.packages("ggtext")
library(ggtext)

save_pdf <- function(filename, plot, width, height, units = "in") {
  showtext_opts(dpi = 72)
  showtext_auto()
  ggsave(filename, plot, width = width, height = height,
         units = units, bg = "white")
  showtext_auto(FALSE)
}

# ---- DATA PREPARATION ----
nat <- df[df$condition == "natural", ]
nat <- nat[order(nat$session), ]
oc  <- df[df$condition == "OC", ]
oc  <- oc[order(oc$session), ]

nat$cycle_day <- ifelse(nat$session >= 11, nat$session - 10, nat$session + 20)
oc$cycle_day  <- oc$session - 30

nat$phase <- cut(nat$cycle_day,
                 breaks = c(0, 4, 11, 14, 30),
                 labels = c("Menses", "Follicular", "Periovulatory", "Luteal"))
oc$phase  <- ifelse(oc$cycle_day >= 14 & oc$cycle_day <= 20, "HFI", "Active")

nat <- nat[order(nat$cycle_day), ]
oc  <- oc[order(oc$cycle_day), ]

nat$log_E2  <- log(nat$`Estradiol (LCMS; pg/mL)`   + 1)
nat$log_P4  <- log(nat$`Progesterone (LCMS; ng/mL)` + 1)
nat$log_LH  <- log(nat$`LH (mIU/ml)`  + 1)
nat$log_FSH <- log(nat$`FSH (mIU/ml)` + 1)

oc$log_E2  <- log(oc$`Estradiol (LCMS; pg/mL)`   + 1)
oc$log_P4  <- log(oc$`Progesterone (LCMS; ng/mL)` + 1)
oc$log_LH  <- log(oc$`LH (mIU/ml)`  + 1)
oc$log_FSH <- log(oc$`FSH (mIU/ml)` + 1)

zscale <- function(x) as.numeric(scale(x))

# Pooled z-scoring (visualisation only)
zscale_pooled <- function(x1, x2) {
  pooled <- c(x1, x2)
  m <- mean(pooled, na.rm=TRUE)
  s <- sd(pooled, na.rm=TRUE)
  list(z1 = (x1 - m) / s, z2 = (x2 - m) / s)
}

pvs_zp  <- zscale_pooled(nat$total_pvs_volume, oc$total_pvs_volume)
psqi_zp <- zscale_pooled(nat$psqi,             oc$psqi)
E2_zp   <- zscale_pooled(nat$log_E2,           oc$log_E2)
P4_zp   <- zscale_pooled(nat$log_P4,           oc$log_P4)
LH_zp   <- zscale_pooled(nat$log_LH,           oc$log_LH)
FSH_zp  <- zscale_pooled(nat$log_FSH,          oc$log_FSH)

nat_plot <- nat %>% mutate(
  pvs_z  = pvs_zp$z1, psqi_z = psqi_zp$z1,
  E2_z   = E2_zp$z1,  P4_z   = P4_zp$z1,
  LH_z   = LH_zp$z1,  FSH_z  = FSH_zp$z1
)

oc_plot <- oc %>% mutate(
  pvs_z  = pvs_zp$z2, psqi_z = psqi_zp$z2,
  E2_z   = E2_zp$z2,  P4_z   = P4_zp$z2,
  LH_z   = LH_zp$z2,  FSH_z  = FSH_zp$z2
)

# Delta data (methodology matches stats.R)
nat_sorted <- nat[order(nat$session), ]
oc_sorted  <- oc[order(oc$session), ]

nat_delta <- data.frame(
  d_pvs     = diff(nat_sorted$total_pvs_volume),
  d_psqi    = diff(nat_sorted$psqi),
  session   = nat_sorted$session[-1],
  cycle_day = nat_sorted$cycle_day[-1],
  phase     = nat_sorted$phase[-1]
)
oc_delta_full <- data.frame(
  d_pvs     = diff(oc_sorted$total_pvs_volume),
  d_psqi    = diff(oc_sorted$psqi),
  session   = oc_sorted$session[-1],
  cycle_day = oc_sorted$cycle_day[-1],
  phase     = oc_sorted$phase[-1]
)

# Exclude pre-registered outlier + following day
z_psqi_oc       <- scale(oc_sorted$psqi)[, 1]
outlier_session <- oc_sorted$session[which(abs(z_psqi_oc) > 3)]

oc_delta <- oc_delta_full[
  !(oc_delta_full$session %in% outlier_session) &
    !(oc_delta_full$session %in% (outlier_session + 1)),
]

# ---- COLOURS ----
hormone_colors <- c(
  "Estradiol"    = "#3D6B94",
  "FSH"          = "#4F9483",
  "LH"           = "#BF8A3D",
  "Progesterone" = "#8B5E83"
)

# Colour-blind-safe palette
phase_colors_nat <- c("Menses"     = "#D55E00",
                      "Follicular"    = "#0072B2",
                      "Periovulatory" = "#009E73",
                      "Luteal"        = "#CC79A7")
phase_colors_oc  <- c("Active" = "#882255",
                      "HFI"    = "#E69F00")

reg_colors <- c(
  delta = "#4A8585",
  e2    = "#3D6B94",
  p4    = "#8B5E83",
  lh    = "#BF8A3D",
  fsh   = "#4F9483"
)

phases_nat <- data.frame(xmin=c(1,12), xmax=c(4,14),
                         label=c("M","O"))
phases_oc  <- data.frame(xmin=14, xmax=20, label="HFI")

# Base font size
TITLE_BASE_SIZE <- 16
LEGEND_TEXT_BASE_SIZE <- 13 * (TITLE_BASE_SIZE / 14)

theme_nature <- function(base_size = TITLE_BASE_SIZE) {
  theme_foundation(base_size = base_size, base_family = "sans") +
    theme(
      text              = element_text(family="sans"),
      plot.title        = element_text(face="bold", size=base_size+4,
                                       hjust=0, margin=margin(b=2),
                                       color="grey10", lineheight=1.2),
      plot.subtitle     = element_text(size=base_size-2, color="grey50",
                                       hjust=0, margin=margin(b=2)),
      axis.title        = element_text(size=base_size, color="black"),
      axis.text         = element_text(size=base_size-2, color="black"),
      axis.ticks        = element_line(color="black", linewidth=0.25),
      axis.ticks.length = unit(2.5, "pt"),
      axis.line         = element_line(color="black", linewidth=0.35),
      panel.background  = element_rect(fill="white", color=NA),
      panel.border      = element_blank(),
      panel.grid        = element_blank(),
      legend.position   = "right",
      legend.key.size   = unit(12, "pt"),
      legend.key.width  = unit(28, "pt"),
      legend.text       = element_text(size=LEGEND_TEXT_BASE_SIZE),
      legend.title      = element_text(size=TITLE_BASE_SIZE, face="bold"),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.spacing.y  = unit(1, "pt"),
      plot.background   = element_rect(fill="white", color=NA),
      plot.margin       = margin(8, 14, 6, 6)
    )
}

# Scale text across figures
panel_w_fig1 <- 18   / (1 + 1 + 0.22)
panel_h_fig1 <- 13.0 / (1 + 1 + 0.12)
panel_w_fig2 <- 18   / 2
panel_h_fig2 <- 7.0  / (1 + 0.08)

ref_panel_area <- panel_w_fig1 * panel_h_fig1

text_scale_fig1 <- 1.15
text_scale_fig2 <- sqrt((panel_w_fig2 * panel_h_fig2) / ref_panel_area)

# Stat-box size scales with panel title size
STAT_SIZE_FRACTION_OF_TITLE <- 0.90
stat_size_for <- function(text_scale) {
  title_pt <- TITLE_BASE_SIZE * text_scale + 4
  (title_pt * STAT_SIZE_FRACTION_OF_TITLE) / 2.845276
}
BASE_POINT_SIZE   <- 3.5
BASE_POINT_STROKE <- 0.5

# Shared PSQI/PVSV axis scale (ref: Figure 2)
REF_PVSV_EXPAND <- 0.18
REF_PSQI_EXPAND <- 0.20

ref_psqi_rng   <- range(nat$psqi, na.rm=TRUE)
ref_psqi_span  <- diff(ref_psqi_rng) * (1 + 0.05 + REF_PSQI_EXPAND)
PSQI_UNITS_PER_INCH <- ref_psqi_span / panel_h_fig2

pvs_lim <- c(floor(min(nat$total_pvs_volume)/50)*50,
             ceiling(max(nat$total_pvs_volume)/50)*50)
ref_pvs_span <- diff(pvs_lim) * (1 + REF_PVSV_EXPAND)
PVSV_UNITS_PER_INCH <- ref_pvs_span / panel_h_fig2

# ---- STATISTICS HELPERS (mirror stats.R) ----
fit_ar <- function(formula_str, data) {
  m0 <- lm(as.formula(formula_str), data = data)
  dw <- dwtest(m0, alternative = "two.sided")
  if (dw$p.value < 0.05) {
    m <- gls(as.formula(formula_str), data = data,
             correlation = corAR1(form = ~time), method = "ML")
    list(model = m, type = "GLS AR(1)")
  } else {
    list(model = m0, type = "OLS")
  }
}

get_coef <- function(fit, term) {
  m <- fit$model
  if (fit$type == "GLS AR(1)") {
    s   <- summary(m)$tTable
    b   <- s[term, "Value"]; se <- s[term, "Std.Error"]
    p   <- s[term, "p-value"]; df_ <- s[term, "DF"]
    ci  <- b + c(-1,1) * qt(0.975, df_) * se
  } else {
    s  <- summary(m)$coefficients
    b  <- s[term, "Estimate"]; se <- s[term, "Std. Error"]
    p  <- s[term, "Pr(>|t|)"]; ci <- confint(m)[term, ]
  }
  list(beta = b, se = se, p = p, ci_lo = ci[1], ci_hi = ci[2])
}

block_permute <- function(x, block_size) {
  n <- length(x)
  x_ext <- c(x, x)
  starts <- sample(seq_len(n), size = ceiling(n / block_size), replace = TRUE)
  blocks <- lapply(starts, function(s) x_ext[s:(s + block_size - 1)])
  unlist(blocks)[seq_len(n)]
}

run_delta_cor <- function(delta_df, block_size = 3, n_perm = 10000) {
  set.seed(42)
  x <- delta_df$d_psqi; y <- delta_df$d_pvs
  ct    <- cor.test(x, y, method = "spearman", exact = FALSE)
  r_obs <- as.numeric(ct$estimate)
  r_perm <- replicate(n_perm,
                      cor(block_permute(x, block_size), y,
                          method = "spearman", use = "complete.obs"))
  p_perm <- (sum(abs(r_perm) >= abs(r_obs)) + 1) / (n_perm + 1)
  list(r = r_obs, p_perm = p_perm)
}

block_bootstrap_ci <- function(delta_df, block_size = 3, n_boot = 10000, conf = 0.95) {
  set.seed(42)
  n <- nrow(delta_df); x <- delta_df$d_psqi; y <- delta_df$d_pvs
  block_id <- ceiling(seq_len(n) / block_size)
  n_blocks <- max(block_id)
  r_boot <- replicate(n_boot, {
    boot_blocks <- sample(seq_len(n_blocks), size = n_blocks, replace = TRUE)
    boot_idx    <- unlist(lapply(boot_blocks, function(b) which(block_id == b)))
    boot_idx    <- boot_idx[seq_len(n)]
    if (length(unique(x[boot_idx])) < 2 || length(unique(y[boot_idx])) < 2) return(NA)
    cor(x[boot_idx], y[boot_idx], method = "spearman")
  })
  r_boot <- r_boot[!is.na(r_boot)]
  c(lo = as.numeric(quantile(r_boot, (1-conf)/2)),
    hi = as.numeric(quantile(r_boot, 1-(1-conf)/2)))
}

sig_star <- function(p) ifelse(p<0.001,"***",ifelse(p<0.01,"**",
                                                    ifelse(p<0.05,"*",ifelse(p<0.10,"\u2020",""))))
# Escaped variant for richtext annotations
sig_star_md <- function(p) ifelse(p<0.001,"\\*\\*\\*",ifelse(p<0.01,"\\*\\*",
                                                             ifelse(p<0.05,"\\*",ifelse(p<0.10,"\u2020",""))))

# Number formatting: mid-dot, 2 sig. digits for p-values
mid_dot <- function(x) gsub("\\.", "&#183;", x)
fmt2 <- function(x) mid_dot(sprintf("%.2f", x))
fmt_p2sig <- function(p) {
  if (p < 0.0001) return("<0\u00b70001")
  rounded <- signif(p, 2)
  if (rounded <= 0) return("<0\u00b70001")
  decimals <- max(0, -floor(log10(rounded)) + 1)
  mid_dot(sprintf(paste0("%.", decimals, "f"), rounded))
}

mid_dot_raw   <- function(x) gsub("\\.", "\u00b7", x)
fmt2_raw      <- function(x) mid_dot_raw(sprintf("%.2f", x))
sig_star_raw  <- function(p) ifelse(p<0.001,"***",ifelse(p<0.01,"**",
                                                         ifelse(p<0.05,"*",ifelse(p<0.10,"\u2020",""))))
fmt_p2sig_raw <- function(p) {
  if (p < 0.0001) return("<0\u00b70001")
  rounded <- signif(p, 2)
  if (rounded <= 0) return("<0\u00b70001")
  decimals <- max(0, -floor(log10(rounded)) + 1)
  mid_dot_raw(sprintf(paste0("%.", decimals, "f"), rounded))
}

# Stat-box text built from grid::textGrob (avoids a gridtext spacing bug)
build_rich_label_grob <- function(segments, x=unit(0.98,"npc"), y=unit(0.95,"npc"),
                                  hjust=c("right","left"),
                                  size=14, sub_scale=0.70, sub_dy_pt=3,
                                  fontface="bold", color="grey20",
                                  fill=scales::alpha("white", 0.75), pad_pt=4) {
  hjust <- match.arg(hjust)
  gps <- lapply(segments, function(seg) {
    fs <- if (isTRUE(seg$sub)) size * sub_scale else size
    grid::gpar(fontsize=fs, fontface=fontface, col=color)
  })
  measure_grobs <- Map(function(seg, gp) grid::textGrob(seg$text, gp=gp), segments, gps)
  widths        <- lapply(measure_grobs, function(g) grid::unit(1, "grobwidth", g))
  total_width   <- Reduce(`+`, widths)

  start_x <- if (hjust == "right") x - total_width - grid::unit(pad_pt, "pt") else x + grid::unit(pad_pt, "pt")

  cur_x <- start_x
  text_grobs <- vector("list", length(segments))
  for (i in seq_along(segments)) {
    seg <- segments[[i]]
    this_y <- if (isTRUE(seg$sub)) y - grid::unit(sub_dy_pt, "pt") else y
    text_grobs[[i]] <- grid::textGrob(seg$text, x=cur_x, y=this_y,
                                      hjust=0, vjust=1, gp=gps[[i]])
    cur_x <- cur_x + widths[[i]]
  }

  bg <- grid::rectGrob(
    x = if (hjust=="right") x else x + grid::unit(pad_pt, "pt")/2 + total_width/2,
    y = y - grid::unit(size*0.35, "pt"),
    width  = total_width + grid::unit(2*pad_pt, "pt"),
    height = grid::unit(size*1.35, "pt"),
    just   = if (hjust=="right") c("right","center") else c("center","center"),
    gp = grid::gpar(fill=fill, col=NA)
  )
  grid::grobTree(bg, do.call(grid::grobTree, text_grobs))
}

make_annotation <- function(p_val, beta) {
  p_str <- ifelse(p_val < 0.0001, "p<0\u00b70001", paste0("p=", fmt_p2sig(p_val)))
  paste0(p_str, sig_star(p_val), "\nb=", fmt2(beta))
}

loess_endpoint <- function(x, y, span=0.35, at=29) {
  ok   <- complete.cases(x, y)
  xok  <- x[ok]; yok <- y[ok]
  fit  <- loess(yok ~ xok, span=span)
  pred <- predict(fit, newdata=data.frame(xok=at[1]))
  as.numeric(pred[1])
}

# Trajectory panel
make_traj <- function(data, phases, title, inset_phase_colors=NULL,
                      pvs_y_override=NULL, psqi_y_override=NULL, text_scale=1) {

  y_pvs_raw  <- loess_endpoint(data$cycle_day, data$pvs_z,  span=0.35, at=29)
  y_psqi_raw <- loess_endpoint(data$cycle_day, data$psqi_z, span=0.35, at=29)

  pvs_label_y  <- if (!is.null(pvs_y_override))  pvs_y_override  else y_pvs_raw
  psqi_label_y <- if (!is.null(psqi_y_override)) psqi_y_override else y_psqi_raw

  if (is.null(pvs_y_override) && is.null(psqi_y_override)) {
    min_gap <- 0.35
    if (!is.na(pvs_label_y) && !is.na(psqi_label_y) &&
        abs(pvs_label_y - psqi_label_y) < min_gap) {
      if (pvs_label_y >= psqi_label_y) {
        pvs_label_y  <- pvs_label_y  + min_gap / 2
        psqi_label_y <- psqi_label_y - min_gap / 2
      } else {
        pvs_label_y  <- pvs_label_y  - min_gap / 2
        psqi_label_y <- psqi_label_y + min_gap / 2
      }
    }
  }

  p <- ggplot(data, aes(x=cycle_day))

  if (nrow(phases) > 0) {
    p <- p +
      geom_rect(data=phases,
                aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf),
                fill="grey93", alpha=1, inherit.aes=FALSE) +
      geom_text(data=phases,
                aes(x=(xmin+xmax)/2, y=5.6, label=label),
                size=5.2 * text_scale, color="black", fontface="bold", inherit.aes=FALSE)
  }

  p <- p +
    geom_ribbon(aes(ymin=0, ymax=E2_z,  fill="Estradiol"),    alpha=0.35) +
    geom_ribbon(aes(ymin=0, ymax=FSH_z, fill="FSH"),          alpha=0.35) +
    geom_ribbon(aes(ymin=0, ymax=LH_z,  fill="LH"),           alpha=0.35) +
    geom_ribbon(aes(ymin=0, ymax=P4_z,  fill="Progesterone"), alpha=0.35) +
    scale_fill_manual(values=hormone_colors,
                      name="Hormones",
                      guide=guide_legend(order=1)) +
    geom_hline(yintercept=0, linetype="dashed",
               color="grey65", linewidth=0.22) +
    geom_smooth(aes(y=psqi_z), color="#B15C3B",
                method="loess", span=0.35, se=FALSE, linewidth=1.0) +
    geom_smooth(aes(y=pvs_z),  color="black",
                method="loess", span=0.35, se=FALSE, linewidth=1.0) +
    annotate("text", x=30.2, y=pvs_label_y,
             label="PVSV",  hjust=0, size=5.0 * text_scale,
             color="black",   fontface="bold") +
    annotate("text", x=30.2, y=psqi_label_y,
             label="1-dPSQI", hjust=0, size=5.0 * text_scale,
             color="#B15C3B", fontface="bold") +
    scale_x_continuous(breaks=seq(1,30,by=3),
                       limits=c(1, 32),
                       expand=c(0,0)) +
    scale_y_continuous(breaks=seq(-3,6,by=1), limits=c(-3,6)) +
    labs(title=title, x="Cycle day", y="z-score (log-transformed)") +
    theme_nature(base_size = TITLE_BASE_SIZE * text_scale) +
    coord_cartesian(clip="off") +
    theme(plot.margin=margin(8, 36, 6, 6))

  if (!is.null(inset_phase_colors)) {
    n         <- length(inset_phase_colors)
    y_pos     <- seq(5.20, by=-0.60, length.out=n)
    x_box_min <- 21.0
    x_box_max <- 29.5
    x_dot     <- 21.8
    x_lbl     <- 22.6

    p <- p +
      annotation_custom(
        grob=rectGrob(gp=gpar(fill="white", col="grey70",
                              lwd=0.5, alpha=0.92)),
        xmin=x_box_min, xmax=x_box_max,
        ymin=min(y_pos) - 0.40, ymax=5.68) +
      annotate("text", x=x_dot, y=5.53,
               label="Phase", size=2.0, fontface="bold",
               color="grey25", hjust=0) +
      annotate("point", x=rep(x_dot, n), y=y_pos,
               color=inset_phase_colors, size=1.5, shape=16) +
      annotate("text", x=rep(x_lbl, n), y=y_pos,
               label=names(inset_phase_colors),
               size=1.9, color="grey25", hjust=0)
  }

  return(p)
}

# Scatter panel
make_scatter <- function(data, x_var, y_var, x_label, y_label,
                         title, line_color, phase_var, phase_colors,
                         p_val, beta, y_limits=NULL, y_expand=0.25,
                         annotation_override=NULL, ann_vjust=1.4,
                         use_richtext=FALSE, text_scale=1, stat_size=5.2) {

  if (!is.null(annotation_override)) {
    annotation <- annotation_override
  } else if (!is.null(p_val)) {
    annotation <- make_annotation(p_val, beta)
  } else {
    annotation <- NULL
  }

  fit_lm  <- lm(as.formula(paste0("`", y_var, "` ~ `", x_var, "`")), data=data)
  x_seq   <- seq(min(data[[x_var]], na.rm=TRUE),
                 max(data[[x_var]], na.rm=TRUE), length.out=200)
  nd      <- setNames(data.frame(x_seq), x_var)
  ci      <- predict(fit_lm, newdata=nd, interval="confidence")
  reg_df  <- data.frame(x=x_seq, fit=ci[,"fit"],
                        lo=ci[,"lwr"], hi=ci[,"upr"])

  p <- ggplot(data, aes(x=.data[[x_var]], y=.data[[y_var]])) +
    geom_ribbon(data=reg_df, aes(x=x, ymin=lo, ymax=hi),
                fill=line_color, alpha=0.12, color=NA,
                inherit.aes=FALSE) +
    geom_line(data=reg_df, aes(x=x, y=fit),
              color=line_color, linewidth=0.9,
              inherit.aes=FALSE) +
    geom_point(aes(fill=.data[[phase_var]]),
               size=BASE_POINT_SIZE * text_scale, alpha=0.88, shape=21,
               stroke=BASE_POINT_STROKE * text_scale, color="white") +
    scale_fill_manual(values=phase_colors, name=NULL, guide="none") +
    scale_x_continuous(breaks=scales::pretty_breaks(n=4),
                       expand=c(0.02, 0)) +
    labs(title=title, x=x_label, y=y_label) +
    theme_nature(base_size = TITLE_BASE_SIZE * text_scale)

  if (!is.null(annotation)) {
    if (use_richtext) {
      p <- p + annotate("richtext", x=Inf, y=Inf,
                        label=annotation,
                        hjust=1, vjust=ann_vjust,
                        size=stat_size, fontface="bold",
                        color="grey20", lineheight=1.3,
                        fill=scales::alpha("white", 0.75), label.color=NA,
                        label.padding=unit(c(3,3,3,3), "pt"))
    } else {
      p <- p + annotate("text", x=Inf, y=Inf,
                        label=annotation,
                        hjust=1, vjust=ann_vjust,
                        size=stat_size, fontface="bold",
                        color="grey20", lineheight=1.3)
    }
  }

  if (!is.null(y_limits)) {
    expanded_lim <- c(y_limits[1], y_limits[2] + diff(y_limits) * y_expand)
    p <- p + scale_y_continuous(limits=expanded_lim,
                                oob=scales::squish,
                                breaks=scales::pretty_breaks(n=4))
  } else {
    y_data <- data[[y_var]]
    y_rng  <- range(y_data, na.rm=TRUE)
    y_pad  <- diff(y_rng) * y_expand
    p <- p + scale_y_continuous(
      limits=c(y_rng[1] - diff(y_rng)*0.05, y_rng[2] + y_pad),
      breaks=scales::pretty_breaks(n=4))
  }
  return(p)
}

make_phase_legend_horiz <- function(phase_colors, text_scale=1) {
  d <- data.frame(
    phase = factor(names(phase_colors), levels=names(phase_colors)),
    x = seq_along(phase_colors), y = 1
  )
  p <- ggplot(d, aes(x=x, y=y, color=phase)) +
    geom_point(size=5.0 * text_scale, shape=16) +
    scale_color_manual(values=phase_colors, name=NULL) +
    guides(color=guide_legend(
      nrow=1, override.aes=list(size=5.0 * text_scale),
      label.position="right",
      keywidth=unit(10, "pt"), keyheight=unit(10, "pt")
    )) +
    theme_void() +
    theme(
      legend.position="bottom", legend.direction="horizontal",
      legend.text=element_text(size=16 * (TITLE_BASE_SIZE/14) * text_scale, color="grey15"),
      legend.key=element_blank(), legend.spacing.x=unit(3, "pt"),
      legend.box.margin=margin(0,0,0,0),
      plot.background=element_rect(fill="white", color=NA)
    )
  cowplot::get_legend(p)
}


# ---- FIGURE 1: TRAJECTORIES + DELTA COUPLING ----

panel_A <- make_traj(nat_plot, phases_nat,
                     "A   Physiological menstrual cycle",
                     inset_phase_colors=NULL,
                     text_scale=text_scale_fig1)
panel_B <- make_traj(oc_plot, phases_oc,
                     "B   Oral contraceptive pill",
                     inset_phase_colors=NULL,
                     pvs_y_override  = -1.20,
                     psqi_y_override = -0.85,
                     text_scale=text_scale_fig1)

leg_traj <- cowplot::get_legend(
  panel_A +
    guides(color="none", fill=guide_legend(order=1)) +
    theme(legend.position="right",
          legend.key.size=unit(14,"pt"),
          legend.key.width=unit(26,"pt"),
          legend.text=element_text(size=LEGEND_TEXT_BASE_SIZE * text_scale_fig1),
          legend.title=element_text(size=TITLE_BASE_SIZE * text_scale_fig1, face="bold"),
          legend.box.margin=margin(0, 0, 0, 8))
)

# Reuse stats.R values if available in session
if (!exists("res_delta_nat")) res_delta_nat <- run_delta_cor(nat_delta, block_size = 3)
if (!exists("res_delta_oc"))  res_delta_oc  <- run_delta_cor(oc_delta,  block_size = 3)
if (!exists("ci_nat"))        ci_nat        <- block_bootstrap_ci(nat_delta, block_size = 3)
if (!exists("ci_oc"))         ci_oc         <- block_bootstrap_ci(oc_delta,  block_size = 3)

seg_nat_delta <- list(
  list(text=paste0("\u03c1=", fmt2_raw(res_delta_nat$r),
                   ", 95% (CI ", fmt2_raw(unname(ci_nat["lo"])), " to ",
                   fmt2_raw(unname(ci_nat["hi"])), "), p"), sub=FALSE),
  list(text="perm", sub=TRUE),
  list(text=paste0("=", fmt_p2sig_raw(res_delta_nat$p_perm),
                   sig_star_raw(res_delta_nat$p_perm)), sub=FALSE)
)
seg_oc_delta <- list(
  list(text=paste0("\u03c1=", fmt2_raw(res_delta_oc$r),
                   ", 95% (CI ", fmt2_raw(unname(ci_oc["lo"])), " to ",
                   fmt2_raw(unname(ci_oc["hi"])), "), p"), sub=FALSE),
  list(text="perm", sub=TRUE),
  list(text=paste0("=", fmt_p2sig_raw(res_delta_oc$p_perm),
                   sig_star_raw(res_delta_oc$p_perm)), sub=FALSE)
)
stat_size_pt_fig1 <- stat_size_for(text_scale_fig1) * 2.845276

# Pin Panel C/D scale to Figure 2
x_rng_raw <- range(c(nat_delta$d_psqi, oc_delta$d_psqi), na.rm=TRUE)
y_rng_raw <- range(c(nat_delta$d_pvs,  oc_delta$d_pvs),  na.rm=TRUE)
x_span_min <- diff(x_rng_raw) * (1 + 0.06)
y_span_min <- diff(y_rng_raw) * (1 + 0.30)

x_span_target <- PSQI_UNITS_PER_INCH * panel_w_fig1
y_span_target <- PVSV_UNITS_PER_INCH * panel_h_fig1
x_span_cd <- max(x_span_min, x_span_target)
y_span_cd <- max(y_span_min, y_span_target)

x_center_cd <- mean(x_rng_raw)
y_center_cd <- mean(y_rng_raw)
x_lim_cd <- c(x_center_cd - x_span_cd/2, x_center_cd + x_span_cd/2)
y_lim_cd <- c(y_center_cd - y_span_cd/2, y_center_cd + y_span_cd/2)

panel_C <- make_scatter(nat_delta, "d_psqi", "d_pvs",
                        "\u03941-dPSQI",
                        "\u0394PVSV (mm\u00b3)",
                        "C   \u03941-dPSQI \u2192 \u0394PVSV\n     Physiological menstrual cycle",
                        reg_colors["delta"], "phase", phase_colors_nat,
                        p_val=NULL, beta=NULL,
                        annotation_override=NULL, use_richtext=FALSE,
                        stat_size=stat_size_for(text_scale_fig1)) +
  coord_cartesian(xlim=x_lim_cd, ylim=y_lim_cd, expand=FALSE) +
  annotation_custom(build_rich_label_grob(seg_nat_delta, y=unit(0.98,"npc"), size=stat_size_pt_fig1),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

panel_D <- make_scatter(oc_delta, "d_psqi", "d_pvs",
                        "\u03941-dPSQI",
                        "\u0394PVSV (mm\u00b3)",
                        "D   \u03941-dPSQI \u2192 \u0394PVSV\n     Oral contraceptive pill",
                        reg_colors["delta"], "phase", phase_colors_oc,
                        p_val=NULL, beta=NULL,
                        annotation_override=NULL, use_richtext=FALSE,
                        stat_size=stat_size_for(text_scale_fig1)) +
  coord_cartesian(xlim=x_lim_cd, ylim=y_lim_cd, expand=FALSE) +
  annotation_custom(build_rich_label_grob(seg_oc_delta, size=stat_size_pt_fig1),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

leg_spacer <- ggplot() + theme_void() +
  theme(plot.background=element_rect(fill="white", color=NA))

f1_A <- panel_A + theme(legend.position="none")
f1_B <- panel_B + theme(legend.position="none")
f1_C <- panel_C + theme(legend.position="none")
f1_D <- panel_D + theme(legend.position="none")

RW <- c(1, 1, 0.26)

fig1_row1 <- cowplot::plot_grid(f1_A, f1_B, leg_traj, ncol=3, rel_widths=RW)
fig1_row2 <- cowplot::plot_grid(f1_C, f1_D, leg_spacer, ncol=3, rel_widths=RW)

leg_strip_nat <- make_phase_legend_horiz(phase_colors_nat, text_scale=text_scale_fig1)
leg_strip_oc  <- make_phase_legend_horiz(phase_colors_oc,  text_scale=text_scale_fig1)

fig1_row3 <- cowplot::plot_grid(
  leg_strip_nat, leg_strip_oc, leg_spacer,
  ncol=3, rel_widths=RW
)

figure1 <- cowplot::plot_grid(
  fig1_row1, fig1_row2, fig1_row3,
  ncol=1, rel_heights=c(1, 1, 0.12)
)

ggsave("figure1.png", figure1, width=18, height=13.0, dpi=600, bg="white")
save_pdf("figure1.pdf", figure1, width=18, height=13.0)


# ---- FIGURE 2: PVSV ~ ESTRADIOL, PSQI ~ LH ----

hormones <- c("log_E2", "log_P4", "log_LH", "log_FSH")

fits_psqi <- lapply(hormones, function(h) fit_ar(paste0("psqi ~ ", h), nat))
names(fits_psqi) <- hormones
fits_pvs  <- lapply(hormones, function(h) fit_ar(paste0("total_pvs_volume ~ ", h), nat))
names(fits_pvs)  <- hormones

p_psqi_all <- sapply(hormones, function(h) get_coef(fits_psqi[[h]], h)$p)
p_pvs_all  <- sapply(hormones, function(h) get_coef(fits_pvs[[h]],  h)$p)
p_psqi_fdr <- p.adjust(p_psqi_all, method="BH")
p_pvs_fdr  <- p.adjust(p_pvs_all,  method="BH")
names(p_psqi_fdr) <- hormones
names(p_pvs_fdr)  <- hormones

cf_pvs_e2  <- get_coef(fits_pvs[["log_E2"]],  "log_E2")
cf_psqi_lh <- get_coef(fits_psqi[["log_LH"]], "log_LH")

p_raw_e2 <- p_pvs_all["log_E2"];  p_fdr_e2 <- p_pvs_fdr["log_E2"]
p_raw_lh <- p_psqi_all["log_LH"]; p_fdr_lh <- p_psqi_fdr["log_LH"]

ann_e2 <- paste0(
  "b=", fmt2(cf_pvs_e2$beta),
  ", 95% (CI ", fmt2(cf_pvs_e2$ci_lo), " to ", fmt2(cf_pvs_e2$ci_hi), ")",
  ", p=", fmt_p2sig(p_raw_e2), sig_star_md(p_raw_e2),
  ", p<sub>FDR</sub>=", fmt_p2sig(p_fdr_e2), sig_star_md(p_fdr_e2)
)
ann_lh <- paste0(
  "b=", fmt2(cf_psqi_lh$beta),
  ", 95% (CI ", fmt2(cf_psqi_lh$ci_lo), " to ", fmt2(cf_psqi_lh$ci_hi), ")",
  ", p=", fmt_p2sig(p_raw_lh), sig_star_md(p_raw_lh),
  ", p<sub>FDR</sub>=", fmt_p2sig(p_fdr_lh), sig_star_md(p_fdr_lh)
)

f2_A <- make_scatter(nat, "log_LH", "psqi",
                     "LH (log)", "1-dPSQI",
                     "A   1-dPSQI ~ LH\n     Physiological menstrual cycle",
                     reg_colors["lh"], "phase", phase_colors_nat,
                     p_val=NULL, beta=NULL,
                     y_expand=0.20,
                     annotation_override=ann_lh,
                     ann_vjust=1.4, use_richtext=TRUE,
                     text_scale=text_scale_fig2,
                     stat_size=stat_size_for(text_scale_fig2))

f2_B <- make_scatter(nat, "log_E2", "total_pvs_volume",
                     "Estradiol (log)", "PVSV (mm\u00b3)",
                     "B   PVSV ~ Estradiol\n     Physiological menstrual cycle",
                     reg_colors["e2"], "phase", phase_colors_nat,
                     p_val=NULL, beta=NULL,
                     y_limits=pvs_lim, y_expand=0.18,
                     annotation_override=ann_e2,
                     ann_vjust=1.4, use_richtext=TRUE,
                     text_scale=text_scale_fig2,
                     stat_size=stat_size_for(text_scale_fig2))

leg_f2 <- make_phase_legend_horiz(phase_colors_nat, text_scale=text_scale_fig2)

f2_row1 <- cowplot::plot_grid(
  f2_A + theme(legend.position="none"),
  f2_B + theme(legend.position="none"),
  ncol=2
)

figure2 <- cowplot::plot_grid(
  f2_row1, leg_f2,
  ncol=1, rel_heights=c(1, 0.08)
)

ggsave("figure2.png", figure2, width=18, height=7.0, dpi=600, bg="white")
save_pdf("figure2.pdf", figure2, width=18, height=7.0)
