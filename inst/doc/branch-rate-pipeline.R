## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4.5, dpi = 140)
library(ontophylo)

## ----workflow, echo=FALSE, fig.width=8.2, fig.height=9, fig.cap="The integrated pipeline. Rate regimes, Neff and Pc are calculated separately for every map replicate before summarizing uncertainty."----
grid::grid.newpage()
box <- function(x, y, width, height, label, fill = "#E1F1FC") {
  grid::grid.roundrect(x, y, width, height, r = grid::unit(0.025, "snpc"),
                       gp = grid::gpar(fill = fill, col = "#555555"))
  grid::grid.text(label, x, y, gp = grid::gpar(fontsize = 10, fontfamily = "sans"))
}
arrow <- function(x0, y0, x1, y1) {
  grid::grid.lines(c(x0, x1), c(y0, y1),
                   arrow = grid::arrow(length = grid::unit(2, "mm")),
                   gp = grid::gpar(col = "#555555"))
}
box(0.5, 0.955, 0.76, 0.06, "Dated tree + observations for each character", "#F3F3F3")
box(0.5, 0.855, 0.84, 0.075, "fit_models()\nCompare Q models independently for each character")
box(0.5, 0.735, 0.84, 0.085,
     "paramo_stm() / read_stm()\nSample conditional histories; retain maps or branch counts")
box(0.5, 0.600, 0.92, 0.10,
     "branch_rate_pipeline()\nCharacter groups + optional per-character normalization\nama_branches_multi_char(): combined branch-rate samples", "#FFF2D5")
box(0.165, 0.435, 0.30, 0.14,
     "Branch samples + summaries\nsummarize_branch_rates_multi()\nplot_branch_colored_tree()")
box(0.5, 0.435, 0.30, 0.14,
     "Marginal rates\nbranch_rate_thru_time_Multi()\nplot_marginal_rates()")
box(0.835, 0.435, 0.30, 0.14,
     "Predefined clade comparisons\nbranch_rate_contrast()\nReplicate differences + ratios")
box(0.5, 0.265, 0.90, 0.10,
     "classify_branch_rates()\nget_enrichment_mask(): quantile regimes per region and map", "#E5F3E7")
box(0.5, 0.100, 0.94, 0.115,
     "phenomic_metrics()\nphenomic_diversity(): Neff+ and Neff-; phenomic_coordination(): Pc\nTime curves, branch coordination and conditional map intervals", "#E5F3E7")
arrow(0.5, 0.924, 0.5, 0.894)
arrow(0.5, 0.816, 0.5, 0.779)
arrow(0.5, 0.691, 0.5, 0.651)
arrow(0.30, 0.549, 0.165, 0.506)
arrow(0.50, 0.549, 0.50, 0.506)
arrow(0.70, 0.549, 0.835, 0.506)
arrow(0.165, 0.364, 0.3, 0.317)
arrow(0.50, 0.214, 0.50, 0.159)

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

## ----phenomic-wrapper---------------------------------------------------------
anatomy <- list(region_A = "char1", region_B = "char2", region_C = "char3")
phenome <- branch_rate_pipeline(mapped, groups = anatomy,
                                 normalize = "characters", n_points = 25,
                                 phenomic = list(lower_q = 0.05, upper_q = 0.95))
head(phenome$phenomic$summary)
head(phenome$phenomic$regimes$branch_summary)
head(phenome$phenomic$regimes$thresholds)

## ----separate-phenomic-functions----------------------------------------------
regimes <- classify_branch_rates(phenome, lower_q = 0.05, upper_q = 0.95)
dim(regimes$accelerated) # branches x regions x maps
dim(regimes$decelerated)
metrics <- phenomic_metrics(phenome, n_points = 25)
stopifnot(isTRUE(all.equal(metrics, phenome$phenomic)))

# Inspect the rate classification for one map replicate.
one_map <- matrix(0, nrow(tree$edge), length(anatomy),
                    dimnames = list(NULL, names(anatomy)))
for (region in names(anatomy)) {
  one_map[, region] <- phenome$groups[[region]]$rate_samples[, 1]
}
accelerated <- get_enrichment_mask(one_map, 0.05, 0.95, direction = "over")
decelerated <- get_enrichment_mask(one_map, 0.05, 0.95, direction = "under")
accelerated_rates <- get_enrichment_mask(one_map, 0.05, 0.95,
                                         direction = "over", return_type = "masked")

## ----exact-phenomic-examples--------------------------------------------------
joint <- rbind(c(1, 1), c(0, 0), c(0, 0), c(0, 0))
separate <- rbind(c(1, 0), c(0, 1), c(0, 0), c(0, 0))
phenomic_coordination(joint)       # Pc = 0.25
phenomic_coordination(separate)    # Pc = 0.25
phenomic_diversity(joint)          # Neff = 0.25
phenomic_diversity(separate)       # Neff = 1

## ----phenomic-plot, fig.cap="Pc calculated per map and then summarized. The band contains pointwise conditional map quantiles."----
pc <- metrics$summary[metrics$summary$metric == "Pc", ]
plot_marginal_rates(pc, ylab = "Phenomic coordination (Pc)")

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

## ----phenomic-validation-summary, echo=FALSE----------------------------------
pilot_file <- system.file("validation", "phenomic-metrics", "results",
                           "pilot-summary.csv", package = "ontophylo")
pilot <- read.csv(pilot_file)
knitr::kable(pilot[, c("scenario", "mean_classified_fraction", "regime_auc",
                       "classified_precision", "Pc_present")], digits = 3,
              col.names = c("Scenario", "Classified fraction", "Regime AUC",
                              "Precision", "Present Pc"))

## ----phenomic-validation-figure, echo=FALSE, out.width="100%", fig.cap="Four fixed simulation pilots. Dotted black curves classify generating rates, dashed orange curves classify the actual histories, and blue curves/bands summarize reconstructed histories (95% conditional map quantiles). All use the same inclusive quantile rule; constant columns are unclassified. These are distinct targets, and the reconstruction does not recover the generating Neff/Pc magnitudes."----
figure <- system.file("validation", "phenomic-metrics", "results",
                      "phenomic-pilot.png", package = "ontophylo")
knitr::include_graphics(figure)

## ----session-info-------------------------------------------------------------
sessionInfo()

