# Run beside run-validation.R after all datasets have completed.
# Rscript --vanilla summarize-validation.R /path/to/results
library(ggplot2)
source("figure-style.R")
arguments <- commandArgs(trailingOnly = TRUE)
output_dir <- "results"
if (length(arguments) == 1L) output_dir <- arguments[1]
stopifnot(length(arguments) <= 1L)
settings <- readRDS(file.path(output_dir, "settings.rds"))
scenarios <- c("homogeneous", "accelerated")
rows <- vector("list", 2L * settings$n_datasets)
index <- 0L
for (scenario in scenarios) {
  for (dataset in seq_len(settings$n_datasets)) {
    index <- index + 1L
    result <- readRDS(file.path(output_dir, sprintf("%s-%02d.rds", scenario, dataset)))
    rows[[index]] <- result$summary
  }
}
results <- do.call(rbind, rows)
stopifnot(nrow(results) == 2L * settings$n_datasets, !anyNA(results$young_branch_auc),
          !anyNA(results$inferred_ratio))
write.csv(results, file.path(output_dir, "recovery-summary.csv"), row.names = FALSE)

summary <- data.frame(Scenario = c("Homogeneous", "10-fold acceleration"),
                      Datasets = settings$n_datasets, Median_AUC = NA_real_,
                      Median_ratio = NA_real_, Excludes_equal = "", Covers_true_ratio = "")
report <- c("Branch-specific rate reconstruction: fixed-tree screening benchmark", "",
             "100 tips; 20 binary characters per dataset; 50 maps per character.",
             "Root age 150 Ma; focal 76 young branches versus 70 other young branches.",
             "Root states have equal frequencies. No pruning or character filtering.",
             "Baseline q=0.002/Myr; accelerated focal q=0.020/Myr.",
             "Each character reconstructed with a separate homogeneous ER fit.",
             "Intervals are conditional map quantiles, not calibrated confidence intervals.", "")
for (i in seq_along(scenarios)) {
  selected <- results[results$scenario == scenarios[i], ]
  n <- nrow(selected)
  summary$Median_AUC[i] <- median(selected$young_branch_auc)
  summary$Median_ratio[i] <- median(selected$inferred_ratio)
  summary$Excludes_equal[i] <- paste0(sum(selected$interval_excludes_equal), "/", n)
  summary$Covers_true_ratio[i] <- paste0(sum(selected$interval_covers_generating_ratio), "/", n)
  exclusion <- binom.test(sum(selected$interval_excludes_equal), n)
  report <- c(report, summary$Scenario[i],
               sprintf("  Datasets: %d", n),
               sprintf("  Median inferred focal/reference ratio: %.4f (range %.4f-%.4f)",
                         median(selected$inferred_ratio), min(selected$inferred_ratio),
                         max(selected$inferred_ratio)),
               sprintf("  Median young-branch AUC: %.4f (range %.4f-%.4f)",
                         median(selected$young_branch_auc), min(selected$young_branch_auc),
                         max(selected$young_branch_auc)),
               sprintf("  Conditional difference interval excludes zero: %d/%d", 
                         sum(selected$interval_excludes_equal), n),
               sprintf("  Exact binomial 95%% interval for this frequency: %.4f-%.4f",
                         exclusion$conf.int[1], exclusion$conf.int[2]),
               sprintf("  Conditional interval entirely positive: %d/%d",
                         sum(selected$positive_interval), n),
               sprintf("  Conditional ratio interval covers the generating ratio: %d/%d",
                         sum(selected$interval_covers_generating_ratio), n),
               sprintf("  Invariant tip-pattern characters retained: %d", sum(selected$invariant_characters)),
               sprintf("  Median correlation with realized branch rates: %.4f",
                         median(selected$branch_rank_correlation)), "")
}
report <- c(report,
  "Interpretation:",
  "AUC assesses localization among young branches, not calibrated branch-wise significance.",
  "Rate magnitude and coverage must be assessed separately from localization.",
  "The null frequency above assesses this fixed design; it is not a universal error guarantee.",
  "Thirty datasets per condition cannot establish precise tail error rates.",
  "The benchmark does not vary tree topology, taxon sampling, resolution or character dependence.",
  "No thresholds were tuned against the results, and no failed datasets were removed.")
write.csv(summary, file.path(output_dir, "scenario-summary.csv"), row.names = FALSE)
writeLines(report, file.path(output_dir, "validation-report.txt"))

results$Scenario <- factor(results$scenario, levels = scenarios,
                           labels = c("Homogeneous", "Accelerated"))
colours <- c(Homogeneous = "#0072B2", Accelerated = "#E69F00")
targets <- data.frame(Scenario = factor(c("Homogeneous", "Accelerated"),
                                        levels = levels(results$Scenario)), ratio = c(1, 10))
jitter <- position_jitter(width = 0.13, height = 0, seed = 19)
auc_plot <- ggplot(results, aes(Scenario, young_branch_auc, colour = Scenario)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "#888888", linewidth = 0.4) +
  geom_point(position = jitter, size = 1.6, alpha = 0.85) +
  scale_colour_manual(values = colours, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = NULL, y = "Young-branch discrimination (AUC)", tag = "a") + theme_nature()
ratio_plot <- ggplot(results, aes(Scenario, inferred_ratio, colour = Scenario)) +
  geom_point(position = jitter, size = 1.6, alpha = 0.85) +
  geom_point(data = targets, aes(Scenario, ratio), inherit.aes = FALSE,
              shape = 18, size = 3, colour = "black") +
  scale_colour_manual(values = colours, guide = "none") +
  scale_y_log10(breaks = c(0.3, 1, 3, 10), limits = c(0.25, 13)) +
  labs(x = NULL, y = "Focal/reference rate ratio (log scale)", tag = "b") + theme_nature()
stopifnot(all(results$inferred_ratio >= 0.25 & results$inferred_ratio <= 13))

draw_figure <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(auc_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(ratio_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::popViewport()
}
grDevices::cairo_pdf(file.path(output_dir, "rate-recovery.pdf"),
                     width = 183 / 25.4, height = 100 / 25.4, family = "Arial", pointsize = 10)
draw_figure()
grDevices::dev.off()
grDevices::png(file.path(output_dir, "rate-recovery.png"), width = 183, height = 100,
               units = "mm", res = 300, type = "cairo", family = "Arial", pointsize = 10)
draw_figure()
grDevices::dev.off()
writeLines(c("Rate-heterogeneity recovery on the fixed 100-tip tree.",
  "Each coloured point is one independently simulated 20-character dataset.",
  "a: Discrimination between focal and reference young branches; the dashed line is chance AUC=0.5.",
  "b: Inferred focal/reference ratio; black diamonds mark the generating ratios (1 and 10).",
  "Inferred ratios compare duration-weighted mean mapped rates. The y axis in b is logarithmic.",
  "Strong discrimination does not imply accurate recovery of the generating fold increase."),
  file.path(output_dir, "figure-caption.txt"))
cat(report, sep = "\n")
