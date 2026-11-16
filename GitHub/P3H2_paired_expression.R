# One-row, three-panel P3H2 paired boxplot figure; combined PNG/PDF only.
# Main inputs: GSE126209 initial 5 FPKM pairs; GSE218035 13 CPM pairs;
# GSE309091 5 pairs of full-matrix normalized counts (previously computed).
# All plotted P values are two-sided paired t tests on log2(Value + 1).
# The exact paired Wilcoxon results are exported alongside them.
# This is a plotting script; it does not fit a count model or recompute size factors.

options(stringsAsFactors = FALSE)

# ============================== Settings ===============================
PROJECT_DIR <- "."
# Optional command-line override: Rscript P3H2_paired_expression.R "another/folder"
cmd_args <- commandArgs(trailingOnly = TRUE)
if (length(cmd_args) >= 1L && nzchar(cmd_args[[1L]])) PROJECT_DIR <- cmd_args[[1L]]

INPUT_DIR <- file.path(PROJECT_DIR, "P3H2_boxplot_inputs")
OUTPUT_DIR <- file.path(PROJECT_DIR, "P3H2_boxplot_outputs")
DATASETS <- c("GSE126209", "GSE218035", "GSE309091")
EXPECTED_PAIRS <- c(GSE126209 = 5L, GSE218035 = 13L, GSE309091 = 5L)
GROUP_COLORS <- c(Normal = "#4C7EB4", Tumor = "#FF4D4F")
FONT_FAMILY <- "Arial"

FONT_SCALE <- 1.3
BASE_SIZE <- 15
TITLE_SIZE <- 12
STRIP_SIZE <- 12
AXIS_TITLE_SIZE <- 12
AXIS_TEXT_SIZE <- 12
LEGEND_TITLE_SIZE <- 12
LEGEND_TEXT_SIZE <- 12
LEGEND_TEXT_SIZE_MANY <- 8
P_VALUE_SIZE <- 12
# ==============================================================
POINT_SIZE <- 2.6
PAIR_LINE_COLOR <- "#9DA58F"
BOX_WIDTH <- 0.56
JITTER_WIDTH <- 0.14
SHOW_X_LABELS <- FALSE
P_VALUE_METHOD <- "paired_t"  # alternatively "wilcoxon_exact"; never chosen by significance
EXPORT_DPI <- 300
PANEL_WIDTH <- 2.7
PANEL_HEIGHT <- 3.8

Y_TITLE_WIDTH <- 0.40
FIGURE_WIDTH <- PANEL_WIDTH * length(DATASETS) + Y_TITLE_WIDTH
FIGURE_HEIGHT <- PANEL_HEIGHT
Y_AXIS_TITLE <- "Expression"
INSTALL_MISSING_PACKAGES <- FALSE

# Register the requested family for Windows graphics/metric calculations.
if (.Platform$OS.type == "windows") {
  do.call(grDevices::windowsFonts,
          stats::setNames(list(grDevices::windowsFont(FONT_FAMILY)), FONT_FAMILY))
}

if (!dir.exists(INPUT_DIR)) stop("Input folder not found: ", INPUT_DIR)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  if (INSTALL_MISSING_PACKAGES) {
    install.packages("ggplot2", repos = "https://cloud.r-project.org")
  } else {
    stop('Please run install.packages("ggplot2") once, then rerun this script.')
  }
}
suppressPackageStartupMessages(library(ggplot2))

natural_patient_order <- function(x) {
  x <- unique(as.character(x))
  if (all(grepl("^[A-Za-z]*[0-9]+$", x))) {
    x[order(as.numeric(gsub("[^0-9]", "", x)), x)]
  } else sort(x)
}

format_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.0001) return(format(p, scientific = TRUE, digits = 2))
  format(signif(p, 2), scientific = FALSE, trim = TRUE)
}

load_dataset <- function(id) {
  path <- file.path(INPUT_DIR, paste0(id, "_P3H2_paired_input.csv"))
  d <- utils::read.csv(path, check.names = FALSE, colClasses = "character", fileEncoding = "UTF-8")
  required <- c("Dataset", "Gene", "Patient", "Sample", "Group", "Value", "Unit", "Batch")
  if (!all(required %in% names(d))) stop("Missing input columns in ", path)
  if (anyNA(d[required]) || any(!nzchar(d$Patient))) stop("Missing sample identifiers in ", id)
  if (!all(d$Dataset == id) || !all(d$Gene == "P3H2")) stop("Unexpected dataset/gene in ", id)
  if (!setequal(unique(d$Group), c("Normal", "Tumor"))) stop("Unexpected group in ", id)
  d$Value <- as.numeric(d$Value)
  if (any(!is.finite(d$Value)) || any(d$Value < 0)) stop("Invalid expression values in ", id)
  if (length(unique(d$Unit)) != 1L || length(unique(d$Batch)) != 1L) stop("Mixed units/batches in ", id)
  if (anyDuplicated(d[c("Patient", "Group")])) stop("Duplicate patient/group in ", id)
  patient_order <- natural_patient_order(d$Patient)
  count_table <- table(factor(d$Patient, levels = patient_order), factor(d$Group, levels = c("Normal", "Tumor")))
  if (!all(count_table == 1L)) stop("Incomplete pairs in ", id)
  if (length(patient_order) != EXPECTED_PAIRS[[id]]) stop("Unexpected pair count in ", id)
  d$Expression <- log2(d$Value + 1)
  d$Patient <- factor(d$Patient, levels = patient_order)
  d$Group <- factor(d$Group, levels = c("Normal", "Tumor"))
  d$Gene <- factor(d$Gene, levels = "P3H2")
  d <- d[order(d$Patient, d$Group), , drop = FALSE]
  # The same deterministic horizontal offset is used for both samples of a pair.
  set.seed(20261111 + match(id, DATASETS))
  offsets <- stats::setNames(stats::runif(length(patient_order), -JITTER_WIDTH, JITTER_WIDTH), patient_order)
  d$x_group <- as.numeric(d$Group)
  d$x_point <- d$x_group + unname(offsets[as.character(d$Patient)])
  d
}

compute_statistics <- function(d) {
  patient_order <- levels(d$Patient)
  n <- d[d$Group == "Normal", ]; t <- d[d$Group == "Tumor", ]
  n <- n[match(patient_order, as.character(n$Patient)), ]; t <- t[match(patient_order, as.character(t$Patient)), ]
  stopifnot(identical(as.character(n$Patient), as.character(t$Patient)))
  differences <- t$Expression - n$Expression
  # Vector pairing is explicit; no formula-based test can reorder the pairs.
  tt <- stats::t.test(t$Expression, n$Expression, paired = TRUE, alternative = "two.sided")
  exact_allowed <- !any(differences == 0) && !anyDuplicated(abs(differences))
  ww <- stats::wilcox.test(t$Expression, n$Expression, paired = TRUE,
                           alternative = "two.sided", exact = exact_allowed, correct = FALSE)
  if (!exact_allowed && P_VALUE_METHOD == "wilcoxon_exact") stop("Exact Wilcoxon unavailable for these differences.")
  data.frame(Dataset = d$Dataset[1L], Gene = "P3H2", Batch = d$Batch[1L],
             n_pairs = length(patient_order), Input_unit = d$Unit[1L], Transform = "log2(Value + 1)",
             mean_paired_log2_difference = mean(differences),
             CI95_low = unname(tt$conf.int[1L]), CI95_high = unname(tt$conf.int[2L]),
             paired_t_statistic = unname(tt$statistic), paired_t_df = unname(tt$parameter),
             paired_t_P = tt$p.value, wilcoxon_statistic = unname(ww$statistic),
             wilcoxon_P = ww$p.value, wilcoxon_exact = exact_allowed,
             tumor_higher_pairs = sum(differences > 0), tumor_lower_pairs = sum(differences < 0),
             stringsAsFactors = FALSE)
}

make_plot <- function(d, s) {
  id <- d$Dataset[1L]; n_patients <- nlevels(d$Patient)
  hues <- seq(15, 375, length.out = n_patients + 1L)[seq_len(n_patients)]
  patient_colors <- stats::setNames(grDevices::hcl(h = hues, c = 100, l = 65), levels(d$Patient))
  spread <- diff(range(d$Expression))
  if (!is.finite(spread) || spread == 0) spread <- 1
  bracket_y <- max(d$Expression) + 0.105 * spread
  label_y <- bracket_y + 0.055 * spread
  shown_p <- if (P_VALUE_METHOD == "paired_t") s$paired_t_P else s$wilcoxon_P
  plot <- ggplot(d, aes(y = Expression)) +
    geom_boxplot(aes(x = x_group, group = Group, colour = Group), fill = "white",
                 width = BOX_WIDTH, linewidth = 0.62, outlier.shape = NA,
                 show.legend = c(colour = TRUE, fill = FALSE)) +
    geom_line(aes(x = x_point, group = Patient), colour = PAIR_LINE_COLOR,
              linewidth = 0.58, alpha = 0.95, show.legend = FALSE) +
    geom_point(aes(x = x_point, fill = Patient, colour = Group), shape = 21,
               size = POINT_SIZE, stroke = 0.58,
               show.legend = c(fill = TRUE, colour = FALSE)) +
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, linewidth = 0.4) +
    annotate("segment", x = 1, xend = 1, y = bracket_y - 0.025 * spread, yend = bracket_y, linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = bracket_y - 0.025 * spread, yend = bracket_y, linewidth = 0.4) +
    annotate("text", x = 1.5, y = label_y, label = format_p(shown_p), size = P_VALUE_SIZE * FONT_SCALE, size.unit = "pt", family = FONT_FAMILY) +
    facet_wrap(~Gene, ncol = 1) +
    scale_colour_manual(name = "Group", values = GROUP_COLORS, breaks = c("Normal", "Tumor"), drop = FALSE) +
    scale_fill_manual(name = "Patient", values = patient_colors, drop = FALSE) +
    scale_x_continuous(limits = c(0.48, 2.52), breaks = if (SHOW_X_LABELS) c(1, 2) else NULL,
                       labels = if (SHOW_X_LABELS) c("Normal", "Tumor") else NULL, expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.055, 0.13))) +
    labs(title = id, x = NULL, y = NULL) +
    guides(colour = guide_legend(order = 1, override.aes = list(fill = "white", linewidth = 0.6)),
           fill = guide_legend(order = 2, override.aes = list(shape = 21, colour = "black", size = 2.4))) +
    theme_bw(base_size = BASE_SIZE * FONT_SCALE, base_family = FONT_FAMILY) +
    theme(panel.grid = element_blank(), panel.border = element_rect(colour = "#333333", linewidth = 0.5),
          strip.background = element_rect(fill = "#D9D9D9", colour = "#333333", linewidth = 0.5),
          strip.text = element_text(size = STRIP_SIZE * FONT_SCALE, margin = margin(5, 5, 5, 5)),
          axis.title.y = element_blank(),
          axis.text = element_text(size = AXIS_TEXT_SIZE * FONT_SCALE, colour = "black"),
          axis.ticks = element_line(colour = "black", linewidth = 0.4),
          plot.title = element_text(size = TITLE_SIZE * FONT_SCALE, hjust = 0.5, margin = margin(b = 3)),
          legend.position = "right", legend.title = element_text(size = LEGEND_TITLE_SIZE * FONT_SCALE),
          legend.text = element_text(size = (if (n_patients > 8L) LEGEND_TEXT_SIZE_MANY else LEGEND_TEXT_SIZE) * FONT_SCALE),
          legend.key.height = grid::unit(if (n_patients > 8L) 0.35 else 0.48, "cm"),
          legend.key.width = grid::unit(0.45, "cm"), legend.spacing.y = grid::unit(0.2, "cm"),
          legend.margin = margin(0, 0, 0, 2), plot.margin = margin(8, 7, 6, 7))
  plot
}

data_list <- stats::setNames(lapply(DATASETS, load_dataset), DATASETS)
statistics <- do.call(rbind, lapply(data_list, compute_statistics))
# This correction concerns the three plotted P3H2 tests, not genome-wide FDR.
statistics$paired_t_BH_across_3_plots <- stats::p.adjust(statistics$paired_t_P, method = "BH")
statistics$plot_P_method <- P_VALUE_METHOD
statistics$plot_P <- if (P_VALUE_METHOD == "paired_t") statistics$paired_t_P else statistics$wilcoxon_P
utils::write.csv(statistics, file.path(OUTPUT_DIR, "P3H2_paired_statistics.csv"), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(do.call(rbind, data_list), file.path(OUTPUT_DIR, "P3H2_plotted_values.csv"), row.names = FALSE, fileEncoding = "UTF-8")

plots <- stats::setNames(lapply(DATASETS, function(id) make_plot(data_list[[id]], statistics[statistics$Dataset == id, ])), DATASETS)

# Force ggplot to build now: catches dropped observations/invalid aesthetics before export.
for (id in DATASETS) {
  built <- ggplot_build(plots[[id]])
  stopifnot(nrow(built$data[[3L]]) == 2L * EXPECTED_PAIRS[[id]])
}

# Each dataset keeps its own y-axis tick values because input scales differ.
# The shared axis title is drawn once, outside the three panels.
draw_combined <- function() {
  grid::grid.newpage()
  grobs <- lapply(plots, ggplotGrob)
  shared_heights <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$heights))
  for (i in seq_along(grobs)) grobs[[i]]$heights <- shared_heights
  widths <- grid::unit.c(grid::unit(Y_TITLE_WIDTH, "in"),
                         grid::unit(rep(1, length(grobs)), "null"))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1L, length(grobs) + 1L, widths = widths)))
  grid::pushViewport(grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L))
  grid::grid.text(Y_AXIS_TITLE, x = 0.55, y = 0.5, rot = 90,
                  gp = grid::gpar(fontfamily = FONT_FAMILY, fontsize = AXIS_TITLE_SIZE * FONT_SCALE))
  grid::popViewport()
  for (i in seq_along(grobs)) {
    grid::pushViewport(grid::viewport(layout.pos.row = 1L, layout.pos.col = i + 1L))
    grid::grid.draw(grobs[[i]])
    grid::popViewport()
  }
  grid::popViewport()
}

export_device <- function(file, kind, width, height, draw) {
  if (kind == "png") {
    grDevices::png(file, width = width, height = height, units = "in", res = EXPORT_DPI,
                   type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"), bg = "white")
  } else if (kind == "pdf") {
    if (capabilities("cairo")) grDevices::cairo_pdf(file, width = width, height = height, family = FONT_FAMILY)
    else grDevices::pdf(file, width = width, height = height, family = "Helvetica", useDingbats = FALSE)
  } else stop("Unsupported output kind: ", kind)
  tryCatch(draw(), finally = grDevices::dev.off())
}

for (kind in c("png", "pdf")) {
  export_device(file.path(OUTPUT_DIR, paste0("P3H2_three_datasets_paired_boxplots.", kind)),
                kind, FIGURE_WIDTH, FIGURE_HEIGHT, draw_combined)
}
utils::capture.output(sessionInfo(), file = file.path(OUTPUT_DIR, "R_sessionInfo.txt"))
cat("Completed. Output folder:\n", normalizePath(OUTPUT_DIR, winslash = "/"), "\n")
print(statistics[c("Dataset", "n_pairs", "paired_t_P", "wilcoxon_P", "plot_P")], row.names = FALSE)

# Documentation:
# https://ggplot2.tidyverse.org/reference/geom_boxplot.html
# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/t.test.html
# https://search.r-project.org/R/refmans/stats/help/wilcox.test.html
