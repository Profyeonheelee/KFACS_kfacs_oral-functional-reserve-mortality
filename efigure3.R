## =============================================================================
## eFigure 3. Sequential pathway attenuation of the association between low
##            functional dentition (<10 units) and all-cause mortality
## =============================================================================
## ---- 0. Packages ------------------------------------------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)
## A UTF-8 locale is required for the en dash used in the confidence intervals.
## Most desktop R installations already provide one, in which case this block
## has no effect; otherwise an ASCII hyphen is substituted.
if (!isTRUE(l10n_info()$`UTF-8`)) {
  for (lc in c("en_US.UTF-8", "C.UTF-8", "C.utf8", "en_GB.UTF-8")) {
    if (suppressWarnings(Sys.setlocale("LC_CTYPE", lc)) != "") break
  }
}
DASH <- if (isTRUE(l10n_info()$`UTF-8`)) "–" else "-"
## ---- 1. Data (Table 3 of the article) ---------------------------------------
dat <- data.frame(
  model = c("Pathway Model 1", "Pathway Model 2", "Pathway Model 3",
            "Pathway Model 4", "Pathway Model 5", "Pathway Model 6A",
            "Pathway Model 6B", "Pathway Model 6C"),
  domain = c("Demographic",
             "+ Socioeconomic and behavioural",
             "+ Weight loss and nutrition",
             "+ Cognition",
             "+ Diagnosed systemic disease burden",
             "+ Physical frailty and laboratory reserve",
             "+ Physical frailty only",
             "+ Laboratory reserve only"),
  beta = c(0.370, 0.329, 0.299, 0.280, 0.279, 0.253, 0.256, 0.275),
  hr   = c(1.45, 1.39, 1.35, 1.32, 1.32, 1.29, 1.29, 1.32),
  lcl  = c(1.18, 1.13, 1.09, 1.07, 1.07, 1.04, 1.05, 1.06),
  ucl  = c(1.78, 1.71, 1.67, 1.64, 1.64, 1.59, 1.60, 1.63),
  p    = c("<.001", ".002", ".006", ".010", ".010", ".020", ".018", ".012"),
  ## Percent attenuation of the log hazard ratio relative to Pathway Model 1,
  ## as reported in Table 3. These were computed from the full-precision
  ## coefficients and are therefore the authoritative values.
  atten = c(NA, 11.0, 19.2, 24.3, 24.6, 31.7, 30.7, 25.6),
  type = c(rep("sequential", 6), rep("branch", 2)),
  stringsAsFactors = FALSE
)
dat$atten_lab <- ifelse(is.na(dat$atten), "Reference", sprintf("%.1f", dat$atten))
dat$hr_lab    <- sprintf("%.2f (%.2f%s%.2f)", dat$hr, dat$lcl, DASH, dat$ucl)
dat$label     <- paste0(dat$model, "    ", dat$domain)
## Integrity check on the tabulated values. The coefficients above are rounded
## to three decimal places, so attenuation recomputed from them agrees with the
## reported values to within about 0.2 percentage points. The tolerance is set
## accordingly: the check detects transcription errors, not rounding.
beta1 <- dat$beta[dat$model == "Pathway Model 1"]
atten_check <- (beta1 - dat$beta) / beta1 * 100
stopifnot(all(abs(atten_check[-1] - dat$atten[-1]) < 0.25))
## ---- 2. Canvas geometry -----------------------------------------------------
## Row 7 is left empty to separate the cumulative models from the branch models.
dat$row <- c(1:6, 8:9)
TOP     <- 9
dat$y   <- TOP - dat$row                       # Pathway Model 1 at the top
X_LAB   <- 0.005                               # model labels, left-aligned
X_F0    <- 0.375; X_F1 <- 0.645                # forest area
X_HR    <- 0.680                               # HR (95% CI), left-aligned
X_P     <- 0.855                               # P value, right-aligned
X_ATT   <- 1.000                               # attenuation, right-aligned
HR_LO   <- 0.97; HR_HI <- 1.90                 # forest x-range
Y_HDR   <- TOP + 0.20                          # header row
Y_AXIS  <- -0.90                               # axis line
Y_TICK  <- -1.32                               # tick labels
## Map a hazard ratio onto the canvas, on a log scale
xf <- function(hr) {
  X_F0 + (log10(hr) - log10(HR_LO)) /
    (log10(HR_HI) - log10(HR_LO)) * (X_F1 - X_F0)
}
dat$x <- xf(dat$hr); dat$xlo <- xf(dat$lcl); dat$xhi <- xf(dat$ucl)
seq_dat <- dat[dat$type == "sequential", ]
seq_dat <- seq_dat[order(seq_dat$row), ]
brks  <- c(1.0, 1.2, 1.4, 1.6, 1.8)
ticks <- data.frame(x = xf(brks), label = sprintf("%.1f", brks))
## ---- 3. Style ---------------------------------------------------------------
## The figure is monochrome by design. Model type is encoded by shape (filled
## circle versus open square) together with the legend, never by colour alone,
## so the figure remains readable in greyscale print and under any form of
## colour-vision deficiency.
INK      <- "#1A1A1A"   # marks and primary text
INK_SOFT <- "#5A5A5A"   # axis text
RULE     <- "#B0B0B0"   # reference line and axis
BASE     <- 9           # base font size (pt)
FAM      <- "sans"      # resolves to Helvetica or Arial on most systems
SZ       <- BASE / .pt  # geom_text size equivalent of BASE
## ---- 4. Figure --------------------------------------------------------------
p <- ggplot() +
  ## reference line at HR = 1.0
  annotate("segment", x = xf(1), xend = xf(1), y = Y_AXIS, yend = TOP - 0.6,
           linetype = "dashed", colour = RULE, linewidth = 0.4) +
  ## line tracing the change in the pooled estimate across cumulative steps
  geom_path(data = seq_dat, aes(x = x, y = y), group = 1,
            colour = INK, linewidth = 0.35, alpha = 0.55) +
  ## confidence intervals and point estimates
  geom_errorbarh(data = dat, aes(y = y, xmin = xlo, xmax = xhi),
                 height = 0, colour = INK, linewidth = 0.55) +
  geom_point(data = dat, aes(x = x, y = y, shape = type, fill = type),
             size = 2.6, colour = INK, stroke = 0.6) +
  ## row labels and numeric columns
  geom_text(data = dat, aes(x = X_LAB, y = y, label = label),
            hjust = 0, size = SZ, colour = INK, family = FAM) +
  geom_text(data = dat, aes(x = X_HR, y = y, label = hr_lab),
            hjust = 0, size = SZ, colour = INK, family = FAM) +
  geom_text(data = dat, aes(x = X_P, y = y, label = p),
            hjust = 1, size = SZ, colour = INK, family = FAM) +
  geom_text(data = dat, aes(x = X_ATT, y = y, label = atten_lab),
            hjust = 1, size = SZ, colour = INK, family = FAM) +
  ## column headers
  annotate("text", x = X_LAB, y = Y_HDR, label = "Pathway model and domain added",
           hjust = 0, fontface = "bold", size = SZ, colour = INK, family = FAM) +
  annotate("text", x = (X_F0 + X_F1) / 2, y = Y_HDR, label = "Pooled estimate",
           hjust = 0.5, fontface = "bold", size = SZ, colour = INK, family = FAM) +
  annotate("text", x = X_HR, y = Y_HDR, label = "HR (95% CI)",
           hjust = 0, fontface = "bold", size = SZ, colour = INK, family = FAM) +
  annotate("text", x = X_P, y = Y_HDR, label = "P value",
           hjust = 1, fontface = "bold", size = SZ, colour = INK, family = FAM) +
  annotate("text", x = X_ATT, y = Y_HDR, label = "Attenuation, %",
           hjust = 1, fontface = "bold", size = SZ, colour = INK, family = FAM) +
  ## x axis for the forest column
  annotate("segment", x = X_F0, xend = X_F1, y = Y_AXIS, yend = Y_AXIS,
           colour = RULE, linewidth = 0.3) +
  annotate("segment", x = ticks$x, xend = ticks$x,
           y = Y_AXIS, yend = Y_AXIS - 0.13, colour = RULE, linewidth = 0.3) +
  annotate("text", x = ticks$x, y = Y_TICK, label = ticks$label,
           size = SZ - 0.15, colour = INK_SOFT, family = FAM, vjust = 1) +
  annotate("text", x = (X_F0 + X_F1) / 2, y = Y_TICK - 0.55,
           label = "Hazard ratio (log scale)",
           size = SZ - 0.15, colour = INK_SOFT, family = FAM, vjust = 1) +
  scale_shape_manual(values = c(sequential = 21, branch = 22),
                     breaks = c("sequential", "branch"),
                     labels = c("Cumulative pathway model",
                                "Branch model (added to Pathway Model 5)"),
                     name = NULL) +
  scale_fill_manual(values = c(sequential = INK, branch = "white"),
                    breaks = c("sequential", "branch"),
                    labels = c("Cumulative pathway model",
                               "Branch model (added to Pathway Model 5)"),
                    name = NULL) +
  coord_cartesian(xlim = c(0, 1), ylim = c(Y_TICK - 1.0, Y_HDR + 0.4),
                  expand = FALSE, clip = "off") +
  theme_void(base_size = BASE, base_family = FAM) +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.text      = element_text(size = BASE - 0.5, colour = INK, family = FAM),
    legend.key.width = unit(14, "pt"),
    legend.margin    = margin(t = 2, b = 0),
    plot.margin      = margin(10, 10, 6, 10),
    plot.background  = element_rect(fill = "white", colour = NA)
  )
## ---- 5. Export --------------------------------------------------------------
W <- 10.2   # width, inches
H <- 4.4    # height, inches
ggsave("eFigure3_pathway_attenuation.pdf",  p, width = W, height = H,
       device = cairo_pdf)
ggsave("eFigure3_pathway_attenuation.png",  p, width = W, height = H,
       dpi = 600, bg = "white")
ggsave("eFigure3_pathway_attenuation.tiff", p, width = W, height = H,
       dpi = 600, bg = "white", compression = "lzw")
## ---- 6. Values plotted ------------------------------------------------------
print(dat[, c("model", "hr", "lcl", "ucl", "p", "atten_lab")], row.names = FALSE)
