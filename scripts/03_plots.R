#!/usr/bin/env Rscript
## 03_plots.R -------------------------------------------------------------
##
## Figures for the support-recovery study.
##
## Two of the editor's requirements are enforced here rather than left to
## chance:
##
##   * the paper is printed in black and white, so nothing is distinguished by
##     colour alone.  Methods are separated by fill pattern (greyscale) and by
##     point shape, and no caption or label refers to a colour;
##   * every piece of text in a figure is at least 10pt.  `base_size` is set
##     to 11 and no element is scaled below 10.
##
## ggplot2 is optional for the rest of the pipeline; if it is missing this
## script exits with a message rather than an error.
##
## Usage:
##   Rscript scripts/03_plots.R [--in=DIR] [--out=DIR]

local({
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(args) > 0L) dirname(sub("^--file=", "", args[1])) else "."
  source(file.path(here, "_bootstrap.R"))
})
root <- bgi_bootstrap()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  message("ggplot2 is not installed; skipping figures.")
  quit(status = 0)
}
library(ggplot2)

parse_flag <- function(args, name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) default else sub(paste0("^--", name, "="), "", hit[1])
}
cli <- commandArgs(trailingOnly = TRUE)
in_dir <- parse_flag(cli, "in", file.path(root, "results", "summaries"))
out_dir <- parse_flag(cli, "out", file.path(root, "results", "figures"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

by_method <- utils::read.csv(
  file.path(in_dir, "support_recovery_by_method.csv"), stringsAsFactors = FALSE)

## Greyscale levels, ordered so that adjacent bars are distinguishable in
## print.  Nothing depends on hue.
method_levels <- c("bgi_sign", "bgi_ci", "bgi_fdr", "pooled_gi",
                   "pooled_gi_bh", "ols", "ols_bh", "icp")
method_labels <- c("BGI (sign)", "BGI (CI)", "BGI (FDR)", "Pooled GI",
                   "Pooled GI + BH", "OLS", "OLS + BH", "ICP")
greys <- grDevices::grey.colors(length(method_levels), start = 0.15,
                                end = 0.92)
names(greys) <- method_levels

by_method$method <- factor(by_method$method, levels = method_levels,
                           labels = method_labels)
by_method$scenario <- factor(by_method$label)

bw_theme <- theme_bw(base_size = 11) +
  theme(
    text = element_text(size = 11),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    strip.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

#' Bar chart of one recovery metric, with Monte Carlo error bars.
metric_plot <- function(metric, se_col, ylab, hline = NA_real_) {
  d <- by_method
  d$value <- d[[metric]]
  d$se <- if (!is.na(se_col) && se_col %in% names(d)) d[[se_col]] else NA_real_

  p <- ggplot(d, aes(x = method, y = value, fill = method)) +
    geom_col(colour = "black", linewidth = 0.3, width = 0.75) +
    scale_fill_manual(values = stats::setNames(greys, method_labels),
                      guide = "none") +
    facet_wrap(~ scenario, ncol = 2) +
    labs(x = NULL, y = ylab) +
    bw_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

  if (any(!is.na(d$se))) {
    p <- p + geom_errorbar(
      aes(ymin = pmax(value - se, 0), ymax = value + se),
      width = 0.25, linewidth = 0.3, na.rm = TRUE)
  }
  if (!is.na(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed",
                        linewidth = 0.4)
  }
  p
}

figs <- list(
  tpr = list(
    plot = metric_plot("tpr", "tpr_se", "Parent recovery (TPR)", 1),
    file = "support_tpr.pdf",
    height = 6
  ),
  fdp = list(
    plot = metric_plot("fdp", "fdp_se",
                       "False discovery proportion",
                       unique(by_method$q_fdr)[1]),
    file = "support_fdp.pdf",
    height = 6
  ),
  fwer = list(
    plot = metric_plot("fwer", "fwer_se",
                       expression(paste("P(", hat(S), " ⊄ pa(Y))")),
                       unique(by_method$alpha)[1]),
    file = "support_fwer.pdf",
    height = 6
  )
)

for (f in figs) {
  ggsave(file.path(out_dir, f$file), f$plot, width = 7, height = f$height,
         device = grDevices::cairo_pdf)
  message("wrote ", file.path(out_dir, f$file))
}

## ---- Calibration ------------------------------------------------------

cal_file <- file.path(in_dir, "calibration_by_scenario.csv")
if (file.exists(cal_file)) {
  cal <- utils::read.csv(cal_file, stringsAsFactors = FALSE)
  long <- do.call(rbind, lapply(
    c("bgi_pred_coverage", "pooled_gi_pred_coverage", "ols_pred_coverage"),
    function(v) {
      data.frame(scenario = cal$label, method = v, coverage = cal[[v]],
                 stringsAsFactors = FALSE)
    }))
  long$method <- factor(long$method,
                        levels = c("bgi_pred_coverage",
                                   "pooled_gi_pred_coverage",
                                   "ols_pred_coverage"),
                        labels = c("BGI", "Pooled GI", "OLS"))

  p <- ggplot(long, aes(x = scenario, y = coverage, shape = method,
                        group = method)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.4) +
    geom_point(size = 2.6, fill = "white", stroke = 0.7) +
    scale_shape_manual(values = c(21, 24, 4), name = NULL) +
    ylim(0, 1) +
    labs(x = NULL, y = "Target-domain predictive coverage") +
    bw_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10))

  ggsave(file.path(out_dir, "predictive_coverage.pdf"), p,
         width = 7, height = 4.2, device = grDevices::cairo_pdf)
  message("wrote ", file.path(out_dir, "predictive_coverage.pdf"))
}

message("Figures written to ", out_dir)
