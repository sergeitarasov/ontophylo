# Run from this directory after loading/installing the development package:
# Rscript --vanilla run-validation.R
# Optional arguments: output directory, datasets per scenario, maps per character.
# The first two datasets form a runtime pilot; settings and seeds are fixed.
# No tree, character, dataset or threshold is selected on recovery performance.
# This benchmark tests the downstream estimator using homogeneous phytools fits.
# The vignette separately demonstrates the fit_models -> paramo_stm interface.

library(ontophylo)
library(ape)
library(phytools)
source("simulate-binary-branches.R")

arguments <- commandArgs(trailingOnly = TRUE)
output_dir <- "results"
n_datasets <- 30L
n_maps <- 50L
if (length(arguments) >= 1) output_dir <- arguments[1]
if (length(arguments) >= 2) n_datasets <- as.integer(arguments[2])
if (length(arguments) >= 3) n_maps <- as.integer(arguments[3])
stopifnot(length(arguments) <= 3, n_datasets >= 1, n_maps >= 2)
first_dataset <- as.integer(Sys.getenv("ONTO_FIRST_DATASET", "1"))
last_dataset <- as.integer(Sys.getenv("ONTO_LAST_DATASET", as.character(n_datasets)))
stopifnot(first_dataset >= 1, last_dataset <= n_datasets, first_dataset <= last_dataset)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

tree <- read.tree("tree.nwk")
branches <- read.csv("branch-regimes.csv")
stopifnot(Ntip(tree) == 100, is.ultrametric(tree),
          identical(unname(tree$edge), unname(as.matrix(branches[, c("parent", "child")]))),
          max(abs(tree$edge.length - branches$duration)) < 1e-8,
          sum(branches$high_rate) == 76,
          sum(branches$young & !branches$high_rate) == 70)
n_char <- 20L
states <- c("0", "1")
scenarios <- c("homogeneous", "accelerated")
focal <- which(branches$high_rate)
reference <- which(branches$young & !branches$high_rate)
young <- which(branches$young)
fast_young <- branches$high_rate[young]
n_fast <- sum(fast_young)
n_slow <- sum(!fast_young)
settings <- list(n_datasets = n_datasets, n_maps = n_maps, n_char = n_char,
                  baseline_q = 0.002, accelerated_q = 0.020,
                  simulation_seed = 61119000L, mapping_seed = 71119000L,
                  interval = "2.5%-97.5% conditional map quantiles",
                  tree_md5 = unname(tools::md5sum("tree.nwk")),
                  code_md5 = unname(tools::md5sum(c("run-validation.R",
                                                   "simulate-binary-branches.R"))))
settings_file <- file.path(output_dir, "settings.rds")
if (file.exists(settings_file)) {
  stopifnot(identical(readRDS(settings_file), settings))
} else {
  saveRDS(settings, settings_file)
}
capture.output(sessionInfo(), file = file.path(output_dir, "session-info.txt"))

for (scenario_id in seq_along(scenarios)) {
  scenario <- scenarios[scenario_id]
  generating_q <- rep(settings$baseline_q, nrow(tree$edge))
  if (scenario == "accelerated") generating_q[focal] <- settings$accelerated_q
  generating_ratio <- if (scenario == "accelerated") 10 else 1

  for (dataset in seq.int(first_dataset, last_dataset)) {
    result_file <- file.path(output_dir, sprintf("%s-%02d.rds", scenario, dataset))
    if (file.exists(result_file)) {
      cat("Already complete:", scenario, dataset, "\n")
      next
    }
    started <- Sys.time()
    simulation_seed <- settings$simulation_seed + 1000L * scenario_id + dataset
    set.seed(simulation_seed)
    roots <- sample(states, n_char, replace = TRUE, prob = c(0.5, 0.5))
    truth <- simulate_binary_branches(tree, generating_q, roots, n_char)
    smaps <- vector("list", n_char)
    names(smaps) <- colnames(truth$tip_states)
    fit_table <- data.frame(character = names(smaps), q = NA_real_,
                            logLik = NA_real_, plateau_gap = NA_real_,
                            invariant = FALSE)
    for (j in seq_len(n_char)) {
      x <- matrix(0, Ntip(tree), 2, dimnames = list(tree$tip.label, states))
      x[cbind(seq_len(Ntip(tree)), match(truth$tip_states[, j], states))] <- 1
      low <- fitMk(tree, x, model = "ER", pi = "equal", q.init = 0.001, logscale = TRUE)
      high <- fitMk(tree, x, model = "ER", pi = "equal", q.init = 0.03, logscale = TRUE)
      fitted <- low
      if (high$logLik > low$logLik + 1e-8) fitted <- high
      if (abs(high$logLik - low$logLik) <= 1e-8 && high$rates < low$rates) fitted <- high
      q <- as.numeric(fitted$rates)
      gap <- as.numeric(fitted$logLik) + Ntip(tree) * log(2)
      fit_table$q[j] <- q
      fit_table$logLik[j] <- fitted$logLik
      fit_table$plateau_gap[j] <- gap
      fit_table$invariant[j] <- length(unique(truth$tip_states[, j])) == 1L
      # An unidentified infinite-rate plateau must not trigger arbitrarily costly maps.
      # Stop and retain the diagnostic dataset rather than dropping a character.
      if (!is.finite(q) || q < 0 || !is.finite(gap) || gap < 1e-6 ||
          fitted$opt_results$convergence != 0) {
        saveRDS(list(scenario = scenario, dataset = dataset, fit_table = fit_table,
                     tip_states = truth$tip_states, settings = settings),
                 file.path(output_dir, sprintf("failed-%s-%02d.rds", scenario, dataset)))
        stop("Fit failed or reached an unidentified high-rate plateau; see saved diagnostic.")
      }
      Q <- matrix(c(-q, q, q, -q), 2, 2, dimnames = list(states, states))
      set.seed(settings$mapping_seed + 100000L * scenario_id + 100L * dataset + j)
      smaps[[j]] <- make.simmap(tree, x, Q = Q, pi = c(0.5, 0.5),
                                nsim = n_maps, tol = 0, message = FALSE)
      for (m in seq_len(n_maps)) {
        mapped_tips <- getStates(smaps[[j]][[m]], "tips")
        stopifnot(identical(unname(mapped_tips[tree$tip.label]),
                            unname(truth$tip_states[, j])))
      }
    }

    comparison <- list(young_clade = list(focal = focal, reference = reference))
    pipeline <- branch_rate_pipeline(smaps, normalize = "characters",
                                      contrasts = comparison, n_points = 61)
    reconstructed <- pipeline$groups$combined
    inferred_comparison <- pipeline$contrasts$combined$young_clade
    # Independent true-history accounting, before inference from tip observations.
    true_result <- branch_rate_pipeline(truth$maps, normalize = "characters",
                                        contrasts = comparison, n_points = 61)
    true_rates <- rowSums(truth$changes_by_branch) / (n_char * tree$edge.length)
    stopifnot(max(abs(true_result$groups$combined$branch_summary$mean - true_rates)) < 1e-10)
    for (j in seq_len(n_char)) {
      stopifnot(all(count_branch_changes(truth$maps[[j]]) == truth$changes_by_branch[, j]))
    }
    inferred_mean <- reconstructed$branch_summary$mean
    scores <- rank(inferred_mean[young], ties.method = "average")
    auc <- (sum(scores[fast_young]) - n_fast * (n_fast + 1) / 2) / (n_fast * n_slow)
    true_scores <- rank(true_rates[young], ties.method = "average")
    true_auc <- (sum(true_scores[fast_young]) - n_fast * (n_fast + 1) / 2) / (n_fast * n_slow)
    ratio <- inferred_comparison$summary[inferred_comparison$summary$metric == "ratio", ]
    difference <- inferred_comparison$summary[
      inferred_comparison$summary$metric == "difference", ]
    ratio_of_means <- mean(inferred_comparison$samples$focal) /
      mean(inferred_comparison$samples$reference)
    observed_ratio <- true_result$contrasts$combined$young_clade$samples$ratio[1]
    # Marginal generating rates and observed event-rate curves are separate targets.
    generating_curve <- branch_rate_thru_time(tree, generating_q, n_points = 61)
    curve <- reconstructed$marginal_summary
    curve$generating_rate <- generating_curve$rate
    curve$realized_rate <- true_result$groups$combined$marginal_summary$mean
    row <- data.frame(scenario = scenario, dataset = dataset,
                       generating_ratio = generating_ratio, realized_ratio = observed_ratio,
                       inferred_ratio = ratio_of_means, ratio_lower = ratio$lower,
                       ratio_upper = ratio$upper,
                       contrast_lower = difference$lower, contrast_upper = difference$upper,
                       interval_excludes_equal = difference$lower > 0 | difference$upper < 0,
                       positive_interval = difference$lower > 0,
                       interval_covers_generating_ratio = ratio$lower <= generating_ratio &
                         ratio$upper >= generating_ratio,
                       young_branch_auc = auc, true_history_auc = true_auc,
                       branch_rank_correlation = cor(inferred_mean, true_rates, method = "spearman"),
                       invariant_characters = sum(fit_table$invariant),
                       minimum_plateau_gap = min(fit_table$plateau_gap),
                       marginal_rmse = sqrt(mean((curve$mean - curve$generating_rate)^2)),
                       seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
    output <- list(summary = row, fit_table = fit_table, curve = curve,
                    branch_summary = reconstructed$branch_summary,
                    branch_truth = data.frame(edge = seq_len(nrow(tree$edge)),
                                              generating_q = generating_q, realized_rate = true_rates),
                    comparison = inferred_comparison, simulation_seed = simulation_seed)
    temporary_file <- paste0(result_file, ".tmp")
    saveRDS(output, temporary_file)
    stopifnot(file.rename(temporary_file, result_file))
    cat(scenario, dataset, "ratio", round(ratio_of_means, 3), "AUC", round(auc, 3),
         "seconds", round(row$seconds, 1), "\n")
    flush.console()
  }
}

rows <- vector("list", length(scenarios) * n_datasets)
index <- 0L
for (scenario in scenarios) {
  for (dataset in seq_len(n_datasets)) {
    result_file <- file.path(output_dir, sprintf("%s-%02d.rds", scenario, dataset))
    if (!file.exists(result_file)) next
    index <- index + 1L
    result <- readRDS(result_file)
    rows[[index]] <- result$summary
  }
}
summary <- do.call(rbind, rows[seq_len(index)])
write.csv(summary, file.path(output_dir, "recovery-summary.csv"), row.names = FALSE)
cat("Completed", nrow(summary), "of", length(scenarios) * n_datasets, "datasets.\n")
print(aggregate(cbind(inferred_ratio, young_branch_auc, interval_excludes_equal,
                       interval_covers_generating_ratio) ~ scenario,
                 data = summary, FUN = mean))
