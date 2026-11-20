# Run with the development version of ontophylo installed or loaded.
# Rscript asr-smoke.R /path/to/output
library(ontophylo)

arguments <- commandArgs(trailingOnly = TRUE)
output_dir <- "asr-validation"
if (length(arguments)) output_dir <- arguments[1]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20261119)
tree <- phytools::pbtree(n = 100, scale = 10)
tree <- ape::reorder.phylo(tree, "cladewise")
n_characters <- 20
true_rate <- 0.1
data <- matrix(NA_character_, 100, n_characters,
                 dimnames = list(tree$tip.label, paste0("char_", seq_len(n_characters))))
true_counts <- matrix(0L, nrow(tree$edge), n_characters)

# Exact symmetric binary process: Poisson jumps flip the state at each jump.
# This simulator is independent of both corHMM and phytools stochastic mapping.
for (character in seq_len(n_characters)) {
  states <- integer(100 + tree$Nnode)
  states[101] <- sample(0:1, 1)
  for (branch in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[branch, 1]
    child <- tree$edge[branch, 2]
    changes <- rpois(1, true_rate * tree$edge.length[branch])
    true_counts[branch, character] <- changes
    states[child] <- (states[parent] + changes) %% 2
  }
  data[, character] <- as.character(states[seq_len(100)])
}

fit_time <- system.time({
  fits <- fit_models(tree, data, models = c("ER", "ARD"), root_p = "equal",
                      nstarts = 3, workers = 2, seed = 1101,
                      outdir = file.path(output_dir, "fits"))
})
stopifnot(all(fits$status$status == "ok"))

mapping_time <- system.time({
  mapped <- paramo_stm(fits, n_stm = 100, batch_size = 25, workers = 2,
                        seed = 1102, store = "branch_counts",
                        outdir = file.path(output_dir, "maps"))
})
stopifnot(all(mapped$status$status == "ok"))
rates <- branch_rates(mapped, n_points = 25)
stopifnot(all(is.finite(rates$rate_samples)))

estimates <- data.frame(character = colnames(data), true_rate = true_rate,
                         selected_model = "", q01 = NA_real_, q10 = NA_real_,
                         true_changes = colSums(true_counts),
                         mean_mapped_changes = NA_real_)
for (i in seq_len(n_characters)) {
  one <- fits$characters[[i]]
  estimates$selected_model[i] <- one$selected_model
  estimates$q01[i] <- one$Q_matrix["0", "1"]
  estimates$q10[i] <- one$Q_matrix["1", "0"]
  counts <- read_stm(mapped, characters = colnames(data)[i], type = "branch_counts")[[1]]
  estimates$mean_mapped_changes[i] <- mean(colSums(counts))
}

write.csv(estimates, file.path(output_dir, "parameter-estimates.csv"), row.names = FALSE)
saveRDS(list(tree = tree, data = data, true_counts = true_counts),
          file.path(output_dir, "simulation.rds"))
saveRDS(fits, file.path(output_dir, "fitted-models.rds"))
saveRDS(mapped, file.path(output_dir, "mapping-index.rds"))
saveRDS(rates, file.path(output_dir, "branch-rates.rds"))
timing <- data.frame(stage = c("fit", "map"),
                      elapsed_seconds = c(fit_time["elapsed"], mapping_time["elapsed"]))
write.csv(timing, file.path(output_dir, "timing.csv"), row.names = FALSE)
print(timing)
print(estimates)
cat("Median selected q01:", median(estimates$q01), "\n")
cat("Median selected q10:", median(estimates$q10), "\n")
cat("Output:", normalizePath(output_dir), "\n")
