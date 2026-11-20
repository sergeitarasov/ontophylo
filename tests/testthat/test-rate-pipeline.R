pipeline_test_map <- function(counts, tree = ape::read.tree(text = "(a:2,b:4);")) {
  tree$maps <- vector("list", nrow(tree$edge))
  for (b in seq_along(counts)) {
    n <- counts[b] + 1
    values <- rep(tree$edge.length[b] / n, n)
    names(values) <- rep(c("0", "1"), length.out = n)
    tree$maps[[b]] <- values
  }
  class(tree) <- c("simmap", "phylo")
  tree
}

test_that("group wrappers preserve memberships, sample pairing and scaling", {
  tree <- ape::read.tree(text = "(a:2,b:2);")
  maps <- list(a = list(pipeline_test_map(c(0, 2), tree), pipeline_test_map(c(2, 0), tree)),
               b = list(pipeline_test_map(c(2, 0), tree), pipeline_test_map(c(0, 4), tree)))
  groups <- list(combined = c("a", "b"), first = "a")
  result <- branch_rate_pipeline(maps, groups, time_points = c(2, 1, 0), time_scale = "age")
  expect_equal(result$groups$combined$rate_samples, branch_rates(maps)$rate_samples)
  expect_equal(result$groups$first$rate_samples, branch_rates(maps["a"])$rate_samples)
  expect_equal(result$group_summary$n_characters, c(2L, 1L))
  expect_equal(nrow(result$branch_summary), 4L)
  expect_equal(nrow(result$marginal_summary), 6L)
  expect_equal(result$memberships, groups)
  mean_rates <- branch_rate_pipeline(maps, groups, normalize = "characters",
                                     time_points = c(2, 1, 0), time_scale = "age")
  expect_equal(mean_rates$groups$combined$rate_samples,
               result$groups$combined$rate_samples / 2)
  expect_equal(mean_rates$groups$combined$marginal_samples,
               result$groups$combined$marginal_samples / 2)
  expect_equal(mean_rates$groups$first$rate_samples, result$groups$first$rate_samples)
  columns <- c("mean", "median", "sd", "lower", "upper")
  expect_equal(mean_rates$groups$combined$branch_summary[, columns],
               result$groups$combined$branch_summary[, columns] / 2)
  expect_equal(mean_rates$groups$combined$marginal_summary[, columns],
               result$groups$combined$marginal_summary[, columns] / 2)
  expect_null(branch_rate_pipeline(maps, marginal = FALSE)$marginal_summary)
})

test_that("invalid group definitions fail rather than silently dropping inputs", {
  maps <- list(a = list(pipeline_test_map(c(0, 2))), b = list(pipeline_test_map(c(2, 0))))
  expect_error(branch_rate_pipeline(unname(maps), marginal = FALSE), "character names")
  expect_error(branch_rate_pipeline(maps, list("a")), "unique names")
  expect_error(branch_rate_pipeline(maps, list(head = c("a", "a"))), "unique, known")
  expect_error(branch_rate_pipeline(maps, list(head = "absent")), "unique, known")
  expect_error(branch_rate_pipeline(maps, list(head = character())), "unique, known")
  expect_error(branch_rate_pipeline(maps, list(head = 1)), "unique, known")
  other <- maps
  other$b <- c(other$b, other$b)
  expect_error(branch_rate_pipeline(other, groups = list(one = "a")), "same number")
  expect_error(branch_rate_pipeline(maps, contrasts = list(x = list(focal = 1))),
               "focal and reference")
})

test_that("clade comparisons preserve exposure weighting and joint map draws", {
  tree <- ape::read.tree(text = "((a:1,b:1):3,c:4);")
  # Two focal branches of different duration: weighted average differs from mean.
  result <- list(tree = tree, rate_samples = cbind(c(0, 4, 1, 1), c(4, 0, 2, 2)))
  comparison <- branch_rate_contrast(result, focal = 1:2, reference = 3:4)
  expect_equal(comparison$samples$focal, c(1, 3))
  expect_equal(comparison$samples$reference, c(1, 2))
  expect_equal(comparison$samples$difference, c(0, 1))
  expect_equal(comparison$samples$ratio, c(1, 1.5))
  expect_equal(comparison$fraction_focal_greater, 0.5)
  expect_equal(comparison$focal_lineage_time, 4)
  expect_equal(comparison$reference_lineage_time, 5)
  equal <- branch_rate_contrast(result, 1:2, 3:4, weighting = "equal")
  expect_equal(equal$samples$focal, c(2, 2))
  expect_error(branch_rate_contrast(result, 1:2, 2:3), "disjoint")
  expect_error(branch_rate_contrast(result, c(1, 1), 2), "unique")
  expect_error(branch_rate_contrast(result, 0, 2), "valid edge")
  expect_error(branch_rate_contrast(result, 1.5, 2), "valid edge")
})

test_that("zero reference rates remain explicit and differences remain available", {
  tree <- ape::read.tree(text = "(a:2,b:2);")
  result <- list(tree = tree, rate_samples = cbind(c(1, 0), c(2, 1)))
  comparison <- branch_rate_contrast(result, 1, 2)
  expect_equal(comparison$samples$ratio, c(NA, 2))
  row <- comparison$summary[comparison$summary$metric == "ratio", ]
  expect_equal(row$n_undefined, 1)
  expect_true(is.na(row$mean))
  expect_true(is.na(row$lower))
  expect_equal(comparison$samples$difference, c(1, 1))
})

test_that("pipeline contrasts match the separate interface", {
  maps <- list(a = list(pipeline_test_map(c(4, 2)), pipeline_test_map(c(6, 4))))
  contrast <- list(focal_clade = list(focal = 1, reference = 2))
  result <- branch_rate_pipeline(maps, contrasts = contrast, marginal = FALSE)
  direct <- branch_rate_contrast(result$groups$combined, 1, 2)
  expect_equal(result$contrasts$combined$focal_clade, direct)
  normalized <- branch_rate_pipeline(maps, contrasts = contrast,
                                      normalize = "characters", marginal = FALSE)
  expect_equal(normalized$contrasts$combined$focal_clade$samples$ratio,
               direct$samples$ratio)
})

test_that("plots retain rates and handle constant colours", {
  tree <- ape::read.tree(text = "(a:2,b:2);")
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  plotted <- plot_branch_colored_tree(tree, c(1, 1), legend = FALSE)
  grDevices::dev.off()
  expect_equal(plotted$value, c(1, 1))
  expect_length(unique(plotted$colour), 1)
  unlink(file)
  summary <- branch_rate_thru_time_Multi(tree, list(c(1, 2), c(2, 1)), n_points = 5)
  plot <- plot_marginal_rates(summary)
  expect_s3_class(plot, "ggplot")
  expect_equal(plot$data, summary)
  expect_error(plot_marginal_rates(rbind(summary, summary)), "unique times")
})

test_that("upstream fits and maps feed the pipeline in memory and from disk", {
  tree <- ape::read.tree(text =
    "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
  data <- cbind(a = c("0", "0", "0", "0", "1", "1", "1", "1"),
                 b = c("0", "0", "1", "1", "0", "0", "1", "1"))
  rownames(data) <- tree$tip.label
  fitted <- fit_models(tree, data, models = "ER", root_p = "equal", nstarts = 0,
                       seed = 61, on_error = "stop")
  mapped <- paramo_stm(fitted, n_stm = 4, batch_size = 2, seed = 62,
                       store = "both", on_error = "stop")
  groups <- list(combined = c("a", "b"), first = "a")
  memory <- branch_rate_pipeline(read_stm(mapped), groups, normalize = "characters")
  direct <- branch_rate_pipeline(mapped, groups, normalize = "characters")
  expect_equal(direct$branch_summary, memory$branch_summary)
  expect_equal(direct$marginal_summary, memory$marginal_summary)
  expect_equal(direct$input_status, mapped$status)

  directory <- tempfile("counts-")
  disk <- paramo_stm(fitted, n_stm = 4, batch_size = 2, seed = 62,
                     store = "branch_counts", outdir = directory, on_error = "stop")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  disk_result <- branch_rate_pipeline(disk, groups, normalize = "characters")
  expect_equal(disk_result$branch_summary, memory$branch_summary)
  expect_equal(disk_result$marginal_summary, memory$marginal_summary)

  mapped$characters$failed <- list(status = "skipped", n_maps = 0L)
  mapped$status <- rbind(mapped$status, data.frame(character = "failed", status = "skipped",
                                                   n_maps = 0L, reason = "test skipped fit"))
  expect_error(branch_rate_pipeline(mapped), "Incomplete character")
  selected <- branch_rate_pipeline(mapped, groups = list(selected = c("a", "b")))
  expect_equal(nrow(selected$input_status), 3)
  expect_true("skipped" %in% selected$input_status$status)
})

test_that("simulated flips have the binary Mk distribution and inherit states", {
  simulator <- system.file("validation", "branch-rates", "simulate-binary-branches.R",
                           package = "ontophylo")
  expect_true(nzchar(simulator))
  source(simulator, local = TRUE)
  tree <- ape::read.tree(text = "(a:1,b:2);")
  set.seed(220611)
  simulation <- simulate_binary_branches(tree, q = c(0.2, 0.4),
                                         root_state = "0", Nchar = 20000,
                                         keep_history = FALSE)
  expected_counts <- c(0.2, 0.8)
  observed_counts <- rowMeans(simulation$changes_by_branch)
  expect_true(all(abs(observed_counts - expected_counts) <
                    5 * sqrt(expected_counts / 20000)))
  expected_tip1 <- (1 - exp(-2 * expected_counts)) / 2
  observed_tip1 <- rowMeans(simulation$tip_states == "1")
  expect_true(all(abs(observed_tip1 - expected_tip1) <
                    5 * sqrt(expected_tip1 * (1 - expected_tip1) / 20000)))
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  histories <- simulate_binary_branches(tree, rep(0.3, 4), "0", 5)
  for (j in 1:5) {
    expect_equal(count_branch_changes(histories$maps[[j]]),
                 unname(histories$changes_by_branch[, j]))
    states <- phytools::getStates(histories$maps[[j]], "tips")
    expect_equal(unname(states[tree$tip.label]), unname(histories$tip_states[, j]))
  }
})
