# Shared style and helpers for all figures (Communications Biology format):
# Arial 8-9 pt, lowercase bold panel letters (10 pt), no titles, lines
# >= 0.5 pt, colour plus shape/linetype coding, vector PDF with embedded
# fonts plus 600 dpi TIFF.

library(dplyr)
library(ggplot2)
library(cowplot)
library(scales)

# ---- SETTINGS ----
FONT_FAMILY <- "Arial"     # "Helvetica" also accepted by the journal
TXT  <- 8                  # pt, axis text, legend text, annotations
TXT2 <- 9                  # pt, axis titles
LETTER_SIZE <- 10          # pt, panel letters
LW_AXIS <- 0.3             # ggplot linewidth units (~0.6 pt)
LW_THIN <- 0.25            # thinnest lines (~0.5 pt)
LW_MAIN <- 0.6             # fitted / trajectory lines
POINT_SIZE   <- 1.6
POINT_STROKE <- 0.3
ALPHA_SIG <- 0.05

pt_to_size <- function(pt) pt / ggplot2::.pt   # text size for geoms

# Font metrics for the null device used by cowplot during layout
if (!FONT_FAMILY %in% names(grDevices::pdfFonts())) {
  do.call(grDevices::pdfFonts, setNames(list(grDevices::pdfFonts()$Helvetica), FONT_FAMILY))
}

# ---- COLOURS, SHAPES, LINETYPES ----
hormone_colors   <- c("Estradiol" = "#3D6B94", "Progesterone" = "#8B5E83",
                      "FSH" = "#3D7A67", "LH" = "#A6762C")
hormone_linetype <- c("Estradiol" = "solid", "Progesterone" = "13",
                      "FSH" = "solid", "LH" = "22")
traj_linetype    <- c("Estradiol" = "solid", "Progesterone" = "13",
                      "FSH" = "42", "LH" = "22")
psqi_color <- "#B15C3B"
pvs_color  <- "black"

phase_colors <- c("Menses" = "#D55E00", "Follicular" = "#0072B2",
                  "Periovulatory" = "#009E73", "Luteal" = "#CC79A7")
phase_shapes <- c("Menses" = 21, "Follicular" = 22,
                  "Periovulatory" = 24, "Luteal" = 23)
phase_colors_oc <- c("Active" = "#882255", "HFI" = "#E69F00")
phase_shapes_oc <- c("Active" = 21, "HFI" = 25)

phase_bands    <- data.frame(xmin = c(1, 12), xmax = c(4, 14), label = c("M", "O"))
phase_bands_oc <- data.frame(xmin = 14, xmax = 20, label = "HFI")

# ---- THEME ----
theme_fig <- function() {
  theme_classic(base_size = TXT, base_family = FONT_FAMILY) +
    theme(
      axis.title        = element_text(size = TXT2, color = "black"),
      axis.text         = element_text(size = TXT,  color = "black"),
      axis.line         = element_line(color = "black", linewidth = LW_AXIS),
      axis.ticks        = element_line(color = "black", linewidth = LW_AXIS),
      axis.ticks.length = unit(1.5, "pt"),
      legend.title      = element_text(size = TXT, face = "bold", color = "black"),
      legend.text       = element_text(size = TXT, color = "black"),
      legend.key.size   = unit(8, "pt"),
      legend.background = element_blank(),
      plot.title        = element_blank(),
      plot.background   = element_rect(fill = "white", color = NA),
      plot.margin       = margin(12, 6, 4, 4)
    )
}

# ---- DATA PREPARATION ----
add_log_hormones <- function(d) {
  d$log_E2  <- log(d$`Estradiol (LCMS; pg/mL)`    + 1)
  d$log_P4  <- log(d$`Progesterone (LCMS; ng/mL)` + 1)
  d$log_LH  <- log(d$`LH (mIU/ml)`  + 1)
  d$log_FSH <- log(d$`FSH (mIU/ml)` + 1)
  d
}

prep_nat <- function(df) {
  nat <- df[df$condition == "natural", ]
  nat <- nat[order(nat$session), ]
  nat$cycle_day <- ifelse(nat$session >= 11, nat$session - 10, nat$session + 20)
  nat$phase <- cut(nat$cycle_day, breaks = c(0, 4, 11, 14, 30),
                   labels = c("Menses", "Follicular", "Periovulatory", "Luteal"))
  add_log_hormones(nat)
}

prep_oc <- function(df) {
  oc <- df[df$condition == "OC", ]
  oc <- oc[order(oc$session), ]
  oc$cycle_day <- oc$session - 30
  oc$phase <- factor(ifelse(oc$cycle_day >= 14 & oc$cycle_day <= 20, "HFI", "Active"),
                     levels = c("Active", "HFI"))
  add_log_hormones(oc)
}

# ---- FORMATTING AND OUTPUT ----
# p-values with two significant digits and leading zero (e.g. "0.034",
# "0.0095"), "< 0.001" below 0.001
fmt_p <- function(p) {
  if (p < 0.001) return("< 0.001")
  paste0("= ", formatC(signif(p, 2), format = "fg", digits = 2, flag = "#"))
}

# Statistics label (plotmath) with white background and no border
stat_layer <- function(label, x, y, hjust = 0) {
  args <- list(geom = "label", x = x, y = y, label = label, parse = TRUE,
               hjust = hjust, vjust = 1, size = pt_to_size(TXT), family = FONT_FAMILY,
               fill = "white", label.padding = unit(2, "pt"))
  if (utils::packageVersion("ggplot2") >= "4.0.0") args$linewidth <- 0 else args$label.size <- 0
  do.call(annotate, args)
}

# Condition/model tag in the top-right corner of the panel, at a fixed
# distance from the panel edges so that it aligns across panels
tag_layer <- function(label) {
  annotation_custom(grid::textGrob(
    label, x = unit(1, "npc") - unit(3, "pt"), y = unit(1, "npc") - unit(3, "pt"),
    hjust = 1, vjust = 1,
    gp = grid::gpar(fontsize = TXT, fontfamily = FONT_FAMILY, col = "black")))
}

# Shaded phase bands with their label inside the top of the band
band_layers <- function(bands) {
  list(
    geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey90", inherit.aes = FALSE),
    geom_text(data = bands, aes(x = (xmin + xmax) / 2, y = Inf, label = label),
              vjust = 1.6, size = pt_to_size(TXT), fontface = "bold", family = FONT_FAMILY,
              inherit.aes = FALSE)
  )
}

phase_legend <- function(colors, shapes, title, nrow = 1) {
  extract_legend(
    ggplot(data.frame(phase = factor(names(colors), levels = names(colors)), x = 1, y = 1),
           aes(x, y, fill = phase, shape = phase)) +
      geom_point(size = 2, stroke = POINT_STROKE, color = "black") +
      scale_fill_manual(values = colors, name = title) +
      scale_shape_manual(values = shapes, name = title) +
      guides(fill = guide_legend(nrow = nrow), shape = guide_legend(nrow = nrow)) +
      theme_fig() +
      theme(legend.position = "bottom", legend.direction = "horizontal"))
}

extract_legend <- function(p) {
  g   <- ggplotGrob(p)
  idx <- grep("^guide-box", g$layout$name)
  idx <- idx[!vapply(g$grobs[idx], inherits, logical(1), "zeroGrob")]
  g$grobs[[idx[1]]]
}

label_panels <- function(..., labels, ncol, align = "hv", axis = "tblr", rel_widths = 1) {
  cowplot::plot_grid(..., ncol = ncol, labels = labels, label_size = LETTER_SIZE,
                     label_fontfamily = FONT_FAMILY, label_fontface = "bold",
                     align = align, axis = axis, rel_widths = rel_widths)
}

# PDF device with embedded fonts: Quartz on macOS (no XQuartz needed),
# Cairo elsewhere, base pdf() with Helvetica as last resort
cairo_ok <- function() {
  isTRUE(capabilities("cairo")) &&
    !inherits(try({
      f <- tempfile(fileext = ".pdf"); grDevices::cairo_pdf(f); grDevices::dev.off(); unlink(f)
    }, silent = TRUE), "try-error")
}
PDF_BACKEND <- if (Sys.info()[["sysname"]] == "Darwin" && isTRUE(capabilities("aqua"))) {
  "quartz"
} else if (suppressWarnings(cairo_ok())) {
  "cairo"
} else {
  "pdf"
}

pdf_device <- function(filename, width, height, ...) {
  switch(PDF_BACKEND,
         quartz = grDevices::quartz(file = filename, type = "pdf", width = width,
                                    height = height, family = FONT_FAMILY, bg = "white"),
         cairo  = grDevices::cairo_pdf(filename, width = width, height = height,
                                       family = FONT_FAMILY,
                                       symbolfamily = grDevices::cairoSymbolFont(FONT_FAMILY, usePUA = FALSE)),
         pdf    = grDevices::pdf(filename, width = width, height = height,
                                 family = "Helvetica", useDingbats = FALSE))
}

save_fig <- function(fig, name, width_mm, height_mm) {
  ggsave(paste0(name, ".pdf"), fig, width = width_mm, height = height_mm, units = "mm",
         device = pdf_device, bg = "white")
  ggsave(paste0(name, ".tiff"), fig, width = width_mm, height = height_mm, units = "mm",
         dpi = 600, device = "tiff", compression = "lzw", bg = "white")
}
