## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4.5, dpi = 140)
library(ontophylo)

## ----workflow, echo=FALSE, fig.width=8.2, fig.height=8, fig.cap="The full pipeline. Branch and marginal summaries retain map replicates until the final statistical summary. Anatomical groups are processed separately."----
grid::grid.newpage()
box <- function(x, y, width, height, label, fill = "#E1F1FC") {
  grid::grid.roundrect(x, y, width, height, r = grid::unit(0.035, "snpc"),
                       gp = grid::gpar(fill = fill, col = "#555555"))
  # Generic sans works with both HTML rendering and R's default PDF device.
  grid::grid.text(label, x, y, gp = grid::gpar(fontsize = 10, fontfamily = "sans"))
}
arrow <- function(x0, y0, x1, y1) {
  grid::grid.lines(c(x0, x1), c(y0, y1),
                   arrow = grid::arrow(length = grid::unit(2, "mm")),
                   gp = grid::gpar(col = "#555555"))
}
box(0.5, 0.94, 0.64, 0.065, "Dated tree + observations for each character", "#F3F3F3")
box(0.5, 0.83, 0.76, 0.075, "fit_models()\nCompare Q models independently for each character")
box(0.5, 0.71, 0.76, 0.075, "paramo_stm()\nSample histories using the fitted Q and root setting")
box(0.5, 0.59, 0.76, 0.075, "read_stm()\nRead a named list of original character-map samples")
box(0.5, 0.45, 0.90, 0.105,
     "branch_rate_pipeline()\nCharacter groups + optional per-character normalization\nama_branches_multi_char(): combined branch-rate samples", "#FFF2D5")
box(0.165, 0.265, 0.30, 0.13,
     "Branch summaries\nsummarize_branch_rates_multi()\nMean, median, SD, quantiles")
box(0.5, 0.265, 0.30, 0.13,
     "Marginal rates\nbranch_rate_thru_time_Multi()\nReplicate curves + summaries")
box(0.835, 0.265, 0.30, 0.13,
     "Predefined clade comparisons\nbranch_rate_contrast()\nReplicate differences + ratios")
box(0.165, 0.08, 0.30, 0.09, "plot_branch_colored_tree()\nBranch-rate tree", "#F3F3F3")
box(0.5, 0.08, 0.30, 0.09, "plot_marginal_rates()\nCurve + mapping interval", "#F3F3F3")
arrow(0.5, 0.906, 0.5, 0.87)
arrow(0.5, 0.791, 0.5, 0.75)
arrow(0.5, 0.671, 0.5, 0.63)
arrow(0.5, 0.551, 0.5, 0.506)
arrow(0.30, 0.396, 0.165, 0.333)
arrow(0.50, 0.396, 0.50, 0.333)
arrow(0.70, 0.396, 0.835, 0.333)
arrow(0.165, 0.198, 0.165, 0.128)
arrow(0.50, 0.198, 0.50, 0.128)

## ----fit-models---------------------------------------------------------------
tree <- ape::read.tree(text =
  "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
characters <- cbind(
  char1 = c("0", "0", "0", "0", "1", "1", "1", "1"),
  char2 = c("0", "0", "1", "1", "0", "0", "1", "1"),
  char3 = c("1", "1", "0", "0", "0", "0", "1", "1"))
rownames(characters) <- tree$tip.label

fits <- fit_models(tree, characters, models = c("ER", "SYM", "ARD"),
                    root_p = "equal", criterion = "AICc", nstarts = 0,
                    seed = 101, on_error = "stop")
fits$status
fits$characters$char1$model_table

mapped <- paramo_stm(fits, n_stm = 25, batch_size = 10,
                      seed = 102, on_error = "stop")
mapped$status
maps <- read_stm(mapped)
names(maps)

## ----high-order-wrapper-------------------------------------------------------
groups <- list(combined = names(maps), region_A = c("char1", "char2"),
                region_B = c("char2", "char3"))
result <- branch_rate_pipeline(maps, groups = groups, normalize = "characters",
                                time_points = seq(3, 0, by = -0.1),
                                time_scale = "age")
result$group_summary
head(result$branch_summary)
head(result$marginal_summary)

combined <- result$groups$combined
dim(combined$rate_samples)       # branches x map replicates
dim(combined$marginal_samples)   # ages x map replicates

## ----direct-mapping-result----------------------------------------------------
direct <- branch_rate_pipeline(mapped, groups, normalize = "characters",
                                time_points = seq(3, 0, by = -0.1),
                                time_scale = "age")
stopifnot(isTRUE(all.equal(direct$branch_summary, result$branch_summary)))

## ----count-and-amalgamate-----------------------------------------------------
counts <- count_branch_changes(maps$char1[[1]])
head(counts)

one_character <- get_branch_rate_sample_smaps(maps$char1)
head(one_character[[1]])

combined_samples <- ama_branches_multi_char(maps)
head(combined_samples[[1]])

single_mean <- summarize_branch_rates(maps$char1, stat = "mean")
combined_mean <- summarize_branch_rates_multi(maps, stat = "mean")
combined_bounds <- summarize_branch_rates_multi(maps, stat = "quantile",
                                                probs = c(0.025, 0.975))
head(single_mean)
head(combined_mean)
combined_bounds[, 1:4, drop = FALSE]

## ----marginal-functions-------------------------------------------------------
one_curve <- branch_rate_thru_time(tree, combined_samples[[1]],
                                   time_points = c(3, 2.5, 2, 1, 0),
                                   time_scale = "age")
one_curve

curves <- branch_rate_thru_time_Multi(tree, combined_samples,
                                      time_points = seq(3, 0, by = -0.1),
                                      time_scale = "age", return_samples = TRUE)
head(curves$summary)

one_combined_character <- branch_rates(maps, n_points = 31)
head(one_combined_character$branch_summary)

## ----clade-contrasts----------------------------------------------------------
focal_tips <- match(c("a", "b", "c", "d"), tree$tip.label)
reference_tips <- match(c("e", "f", "g", "h"), tree$tip.label)
focal <- which(tree$edge[, 2] %in% focal_tips)
reference <- which(tree$edge[, 2] %in% reference_tips)
comparison <- branch_rate_contrast(combined, focal, reference)
comparison$summary
comparison$fraction_focal_greater

comparisons <- list(left_vs_right = list(focal = focal, reference = reference))
with_contrasts <- branch_rate_pipeline(maps, groups,
                                       normalize = "characters", contrasts = comparisons)
with_contrasts$contrasts$combined$left_vs_right$summary

## ----branch-plot, fig.cap="Combined-character branch rates in the example, expressed per character per tree time unit."----
plot_branch_colored_tree(combined$tree, combined$branch_summary$mean,
                         legend_title = "Rate per character", edge.width = 3,
                         cex = 0.8)

## ----marginal-plot, fig.cap="Marginal rates over sampled lineages. Shading is the pointwise central 95% interval across map replicates, conditional on the fitted models and tree."----
plot_marginal_rates(combined$marginal_summary,
                     ylab = "Mapped changes per character per time unit")

## ----validation-summary, echo=FALSE-------------------------------------------
summary_file <- system.file("validation", "branch-rates", "results",
                            "scenario-summary.csv", package = "ontophylo")
stopifnot(nzchar(summary_file))
validation <- read.csv(summary_file)
knitr::kable(validation, digits = 3,
              col.names = c("Condition", "Datasets", "Median AUC", "Median ratio",
                              "Interval excludes equal rates", "Interval covers true ratio"))

## ----validation-figure, echo=FALSE, out.width="100%", fig.cap="Validation on the fixed 100-tip tree. Each point represents one independently simulated character dataset. AUC describes discrimination between focal and reference young branches; it is not a branch-wise significance test."----
figure <- system.file("validation", "branch-rates", "results",
                      "rate-recovery.png", package = "ontophylo")
stopifnot(nzchar(figure))
knitr::include_graphics(figure)

## ----session-info-------------------------------------------------------------
sessionInfo()

