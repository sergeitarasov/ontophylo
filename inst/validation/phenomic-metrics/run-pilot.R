# Run with the new ontophylo installed: Rscript run-pilot.R new-results
# Four prespecified smoke tests, not a power study or a calibrated hypothesis test.
library(ontophylo)
source(system.file('validation', 'branch-rates', 'simulate-binary-branches.R',
                   package = 'ontophylo'))
tree <- ape::read.tree(system.file('validation', 'branch-rates', 'tree.nwk',
                                   package = 'ontophylo'))
branches <- read.csv(system.file('validation', 'branch-rates', 'branch-regimes.csv',
                                 package = 'ontophylo'))
stopifnot(all(tree$edge == as.matrix(branches[, c('parent', 'child')])))
arguments <- commandArgs(trailingOnly = TRUE)
output_dir <- 'results'
if (length(arguments) == 1) output_dir <- arguments[1]
stopifnot(length(arguments) <= 1)
if (file.exists(file.path(output_dir, 'settings.rds'))) {
  stop('Results already exist. Supply a new output directory; existing runs are preserved.')
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
scenarios <- c('homogeneous', 'coordinated', 'mosaic', 'decelerated')
groups <- list(region_A = paste0('CH', 1:10), region_B = paste0('CH', 11:20))
settings <- list(n_tips = 100, n_characters = 20, n_maps = 25,
                  baseline_q = 0.002, high_q = 0.020, low_q = 0.0002,
                  lower_q = 0.05, upper_q = 0.95,
                  seed = 20261120L, scenarios = scenarios,
                  model = 'ER', root_p = 'equal', nstarts = 1L)
saveRDS(settings, file.path(output_dir, 'settings.rds'))
rows <- vector('list', length(scenarios))
for (i in seq_along(scenarios)) {
  scenario <- scenarios[i]
  started <- Sys.time()
  q <- matrix(settings$baseline_q, nrow(tree$edge), 2,
               dimnames = list(NULL, names(groups)))
  if (scenario == 'coordinated') q[branches$high_rate, ] <- settings$high_q
  if (scenario == 'mosaic') {
    q[branches$high_rate, 1] <- settings$high_q
    q[branches$young & !branches$high_rate, 2] <- settings$high_q
  }
  if (scenario == 'decelerated') q[branches$high_rate, ] <- settings$low_q
  truth_maps <- vector('list', 20)
  names(truth_maps) <- paste0('CH', 1:20)
  observations <- matrix('', 100, 20, dimnames = list(tree$tip.label, names(truth_maps)))
  for (j in 1:2) {
    set.seed(settings$seed + i * 100L + j)
    roots <- sample(c('0', '1'), 10, replace = TRUE)
    simulation <- simulate_binary_branches(tree, q[, j], roots, 10)
    ids <- ((j - 1) * 10 + 1):(j * 10)
    observations[, ids] <- simulation$tip_states
    for (k in 1:10) truth_maps[[ids[k]]] <- simulation$maps[[k]]
  }
  fits <- fit_models(tree, observations, models = 'ER', root_p = 'equal',
                      nstarts = 1, seed = settings$seed + i * 1000L,
                      on_error = 'stop', control = list(constant = 'error'))
  saveRDS(fits$status, file.path(output_dir, paste0(scenario, '-fit-status.rds')))
  stopifnot(all(fits$status$status == 'ok'))
  mapped <- paramo_stm(fits, n_stm = settings$n_maps, batch_size = 25,
                        store = 'branch_counts', seed = settings$seed + i * 10000L,
                        on_error = 'stop')
  inferred <- branch_rate_pipeline(mapped, groups, normalize = 'characters',
                                    n_points = 61, phenomic = list())
  realized <- branch_rate_pipeline(truth_maps, groups, normalize = 'characters',
                                    n_points = 61, phenomic = list())
  generating_result <- list(groups = list(), memberships = groups)
  for (j in 1:2) {
    generating_result$groups[[names(groups)[j]]] <- list(tree = tree,
                                                        rate_samples = matrix(q[, j], ncol = 1))
  }
  generating <- phenomic_metrics(generating_result, n_points = 61)
  branch_table <- inferred$phenomic$regimes$branch_summary
  selected_direction <- if (scenario == 'decelerated') 'decelerated_frequency' else 'accelerated_frequency'
  known <- logical(nrow(branch_table))
  for (j in 1:2) {
    ids <- which(branch_table$group == names(groups)[j])
    if (scenario == 'decelerated') known[ids] <- q[, j] < settings$baseline_q
    if (scenario != 'decelerated') known[ids] <- q[, j] > settings$baseline_q
  }
  score <- branch_table[[selected_direction]]
  auc <- precision <- NA_real_
  if (any(known) && any(!known)) {
    ranks <- rank(score, ties.method = 'average')
    n_true <- sum(known)
    n_other <- sum(!known)
    auc <- (sum(ranks[known]) - n_true * (n_true + 1) / 2) / (n_true * n_other)
    if (sum(score) > 0) precision <- sum(score[known]) / sum(score)
  }
  pc <- inferred$phenomic$summary
  pc <- pc[pc$metric == 'Pc', ]
  rows[[i]] <- data.frame(scenario = scenario, generating_regime_fraction = mean(known),
                           mean_classified_fraction = mean(score),
                           regime_auc = auc, classified_precision = precision,
                           Pc_present = pc$mean[pc$time == 0],
                           elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = 'secs')))
  saveRDS(list(settings = settings, observations = observations, generating_q = q,
                fitted = fits, mapped = mapped, inferred = inferred,
                realized = realized, generating = generating),
            file.path(output_dir, paste0(scenario, '.rds')))
  print(rows[[i]])
  write.csv(do.call(rbind, rows[seq_len(i)]), file.path(output_dir, 'pilot-summary.csv'),
             row.names = FALSE)
}
capture.output(sessionInfo(), file = file.path(output_dir, 'session-info.txt'))
