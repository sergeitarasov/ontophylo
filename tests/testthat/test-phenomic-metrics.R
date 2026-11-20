phenomic_test_result <- function() {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  first <- c(0, 10, 0, 0, 0, 0)
  second <- c(0, 0, 10, 0, 0, 0)
  list(groups = list(A = list(tree = tree, rate_samples = cbind(first, first)),
                     B = list(tree = tree, rate_samples = cbind(first, second))),
       memberships = list(A = "char_A", B = "char_B"))
}

test_that("classification uses separate region quantiles and explicit ties", {
  rates <- cbind(A = 0:4, B = c(40, 30, 20, 10, 0))
  rownames(rates) <- letters[1:5]
  high <- get_enrichment_mask(rates, 0.25, 0.75, "over")
  low <- get_enrichment_mask(rates, 0.25, 0.75, "under")
  expect_equal(unname(high), cbind(c(0L, 0L, 0L, 1L, 1L), c(1L, 1L, 0L, 0L, 0L)))
  expect_equal(unname(low), cbind(c(1L, 1L, 0L, 0L, 0L), c(0L, 0L, 0L, 1L, 1L)))
  expect_equal(dimnames(high), dimnames(rates))
  both <- get_enrichment_mask(rates, 0.25, 0.75, "both")
  expect_equal(both, pmax(high, low))
  strict <- get_enrichment_mask(rates, 0.25, 0.75, "over", ties = "exclude")
  expect_equal(unname(strict), cbind(c(0L, 0L, 0L, 0L, 1L), c(1L, 0L, 0L, 0L, 0L)))
  expect_equal(get_enrichment_mask(rates, 0.25, 0.75, "over", "masked"), rates * high)
  scaled <- rates
  scaled[, 1] <- scaled[, 1] * 100
  scaled[, 2] <- scaled[, 2] / 100
  expect_equal(get_enrichment_mask(scaled, 0.25, 0.75, "over"), high)
})

test_that("constant and zero rates cannot silently imply acceleration", {
  rates <- cbind(zero = rep(0, 5), constant = rep(2, 5))
  expect_true(all(get_enrichment_mask(rates, direction = "over") == 0))
  expect_true(all(get_enrichment_mask(rates, direction = "under") == 0))
  expect_true(all(get_enrichment_mask(rates, constant = "include") == 1))
  expect_true(all(get_enrichment_mask(rates, ties = "exclude", constant = "include") == 0))
  expect_equal(get_enrichment_mask(matrix(2, 1, 1)), matrix(0L, 1, 1))
  many_zeros <- matrix(c(rep(0, 19), 1), ncol = 1)
  expect_equal(sum(get_enrichment_mask(many_zeros, 0.05, 0.95, "under")), 19)
  expect_equal(sum(get_enrichment_mask(many_zeros, 0.05, 0.95, "under",
                                      ties = "exclude")), 0)
})

test_that("regime input failures are explicit", {
  rates <- matrix(1:4, 2)
  expect_error(get_enrichment_mask(as.data.frame(rates)), "matrix")
  expect_error(get_enrichment_mask(matrix(numeric(), 0, 1)), "matrix")
  expect_error(get_enrichment_mask(matrix(c(0, NA), 2)), "finite")
  expect_error(get_enrichment_mask(-rates), "nonnegative")
  expect_error(get_enrichment_mask(rates, 0.9, 0.1), "lower_q")
  expect_error(get_enrichment_mask(rates, NA, 0.95), "lower_q")
  expect_error(get_enrichment_mask(rates, 0, 1.1), "lower_q")
  expect_error(get_enrichment_mask(rates, c(0, 0.1), 1), "lower_q")
})

test_that("Neff is whole-pattern diversity rather than region diversity", {
  joint <- rbind(c(1, 1), c(0, 0), c(0, 0), c(0, 0))
  separate <- rbind(c(1, 0), c(0, 1), c(0, 0), c(0, 0))
  expect_equal(phenomic_coordination(joint)$Pc, 0.25)
  expect_equal(phenomic_coordination(separate)$Pc, 0.25)
  expect_equal(phenomic_diversity(joint)$Neff, 0.25)
  expect_equal(phenomic_diversity(separate)$Neff, 1)
  expect_equal(phenomic_diversity(separate)$entropy, log(2))
  expect_equal(phenomic_diversity(separate)$fraction, 0.5)
  expect_equal(phenomic_diversity(separate)$n_patterns, 2)
  skewed <- rbind(c(1, 0), c(1, 0), c(1, 0), c(0, 1), c(0, 0))
  expected_entropy <- -(0.75 * log(0.75) + 0.25 * log(0.25))
  expected_neff <- 0.8 * exp(expected_entropy)
  expect_equal(phenomic_diversity(skewed)$Neff, expected_neff)
  expect_equal(phenomic_diversity(skewed[5:1, 2:1]), phenomic_diversity(skewed))
  expect_equal(phenomic_diversity(skewed > 0), phenomic_diversity(skewed))
  expect_equal(phenomic_coordination(skewed)$branch, c(0.5, 0.5, 0.5, 0.5, 0))
  expect_equal(phenomic_diversity(matrix(1, 1, 1))$Neff, 1)
  expect_equal(phenomic_coordination(matrix(1, 1, 1))$Pc, 1)
  empty <- phenomic_diversity(matrix(0, 3, 2))
  expect_equal(empty$Neff, 0)
  expect_equal(empty$effective_patterns, 0)
  expect_equal(empty$fraction, 0)
  expect_equal(empty$n_lineages, 3)
  expect_equal(phenomic_coordination(matrix(0, 3, 2))$Pc, 0)
  expect_error(phenomic_diversity(matrix(NA, 3, 2)), "0/1 matrix")
  expect_error(phenomic_coordination(matrix(2, 3, 2)), "0/1 matrix")
  expect_error(phenomic_diversity(matrix(0, 0, 2)), "nonempty")
  expect_error(phenomic_coordination(c(0, 1)), "matrix")
})

test_that("regime arrays preserve map pairing and region order", {
  result <- phenomic_test_result()
  # Map-column names must describe matched draws, not source vector names.
  colnames(result$groups$B$rate_samples) <- colnames(result$groups$A$rate_samples)
  regimes <- classify_branch_rates(result)
  expect_equal(dim(regimes$accelerated), c(6L, 2L, 2L))
  expect_equal(regimes$accelerated[2, , 1], c(A = 1L, B = 1L))
  expect_equal(regimes$accelerated[2, , 2], c(A = 1L, B = 0L))
  expect_equal(regimes$accelerated[3, , 2], c(A = 0L, B = 1L))
  expect_equal(regimes$thresholds$upper, rep(7.5, 4))
  expect_equal(regimes$thresholds$n_accelerated, rep(1, 4))
  expect_equal(regimes$thresholds$n_decelerated, rep(5, 4))
  expect_equal(regimes$branch_summary$accelerated_frequency,
               c(0, 1, 0, 0, 0, 0, 0, 0.5, 0.5, 0, 0, 0))
  reverse <- classify_branch_rates(result, groups = c("B", "A"))
  expect_equal(reverse$accelerated, regimes$accelerated[, c("B", "A"), ])
  one <- classify_branch_rates(result, "A")
  expect_equal(dim(one$accelerated), c(6L, 1L, 2L))
})

test_that("nonlinear metrics are calculated before averaging maps", {
  result <- phenomic_test_result()
  colnames(result$groups$B$rate_samples) <- colnames(result$groups$A$rate_samples)
  metrics <- phenomic_metrics(result, time_points = c(0, 2, 1, 0.2), time_scale = "age")
  expect_equal(unname(metrics$samples$Neff_plus[1, ]), c(0.25, 1))
  expect_equal(unname(metrics$samples$Pc[1, ]), c(0.25, 0.25))
  expect_equal(unname(metrics$samples$Neff_minus[1, ]), c(0.75, sqrt(8)))
  plus <- metrics$summary[metrics$summary$metric == "Neff_plus", ]
  expect_equal(plus$time, c(0, 2, 1, 0.2))
  expect_equal(plus$n_branches, c(4L, 2L, 4L, 4L))
  expect_equal(plus$mean, c(0.625, 0, 0.625, 0.625))
  expect_equal(plus$lower[1], 0.26875)
  expect_equal(plus$upper[1], 0.98125)
  expect_equal(unname(metrics$samples$Neff_plus[3, ]), c(0.25, 1))
  elapsed <- phenomic_metrics(result, time_points = c(2, 0, 1, 1.8))
  expect_equal(elapsed$samples, metrics$samples)
  rates <- branch_rate_thru_time(result$groups$A$tree, rep(1, 6),
                                 time_points = c(0, 2, 1, 0.2), time_scale = "age")
  expect_equal(plus$n_branches, rates$n_branches)
  expect_equal(metrics$branch_coordination_summary$mean, c(0, 0.75, 0.25, 0, 0, 0))
  # A mean-rate-first analysis would classify only A on branch 2: Neff = 0.25.
  averaged <- result
  for (g in names(averaged$groups)) {
    averaged$groups[[g]]$rate_samples <- matrix(rowMeans(result$groups[[g]]$rate_samples), 6)
  }
  mean_first <- phenomic_metrics(averaged, time_points = 0, time_scale = "age")
  expect_false(isTRUE(all.equal(mean_first$summary$mean[1], plus$mean[1])))
  expect_true(all(is.na(mean_first$summary$sd)))
})

test_that("constant regions stay in the denominator and overlap is reported", {
  result <- phenomic_test_result()
  colnames(result$groups$B$rate_samples) <- colnames(result$groups$A$rate_samples)
  result$groups$B$rate_samples[,] <- 0
  metrics <- phenomic_metrics(result, time_points = 0, time_scale = "age")
  expect_equal(unname(metrics$samples$Pc), matrix(1/8, 1, 2))
  expect_equal(unname(metrics$samples$Neff_plus), matrix(1/4, 1, 2))
  expect_equal(metrics$regimes$thresholds$constant, c(FALSE, TRUE, FALSE, TRUE))
  legacy <- classify_branch_rates(result, constant = "include")
  expect_true(all(legacy$accelerated[, "B", ] == 1))
  expect_true(all(legacy$decelerated[, "B", ] == 1))
  result$memberships$B <- "char_A"
  expect_warning(regimes <- classify_branch_rates(result), "share characters")
  expect_equal(regimes$shared_characters, "char_A")
})

test_that("mismatched inputs and times cannot produce phenomic summaries", {
  result <- phenomic_test_result()
  expect_error(classify_branch_rates(result), "map columns")
  colnames(result$groups$B$rate_samples) <- colnames(result$groups$A$rate_samples)
  expect_error(classify_branch_rates(list()), "pipeline result")
  expect_error(classify_branch_rates(result, c("A", "A")), "unique")
  expect_error(classify_branch_rates(result, "absent"), "known")
  expect_error(phenomic_metrics(result, time_points = 3), "tree height")
  expect_error(phenomic_metrics(result, n_points = 1), "at least two")
  wrong <- result
  wrong$groups$B$tree$tip.label <- rev(wrong$groups$B$tree$tip.label)
  expect_error(classify_branch_rates(wrong), "same tree")
  wrong <- result
  wrong$groups$B$rate_samples <- wrong$groups$B$rate_samples[, 1, drop = FALSE]
  expect_error(classify_branch_rates(wrong), "map columns")
  for (g in names(result$groups)) result$groups[[g]]$tree$edge.length[2] <- 0.5
  expect_s3_class(classify_branch_rates(result), "ontophylo_regimes")
  expect_error(phenomic_metrics(result), "ultrametric")
})

test_that("the full map pipeline includes phenomic metrics only when requested", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  # Make alternating-state maps directly so expected branch counts are exact.
  make_maps <- function(counts) {
    maps <- vector("list", ncol(counts))
    for (i in seq_len(ncol(counts))) {
      map <- tree
      map$maps <- vector("list", nrow(tree$edge))
      for (b in seq_len(nrow(tree$edge))) {
        n <- counts[b, i] + 1
        segments <- rep(tree$edge.length[b] / n, n)
        names(segments) <- rep(c("0", "1"), length.out = n)
        map$maps[[b]] <- segments
      }
      class(map) <- c("simmap", "phylo")
      maps[[i]] <- map
    }
    maps
  }
  maps <- list(a = make_maps(cbind(c(0, 10, 0, 0, 0, 0), c(0, 10, 0, 0, 0, 0))),
               b = make_maps(cbind(c(0, 10, 0, 0, 0, 0), c(0, 0, 10, 0, 0, 0))))
  groups <- list(combined = c("a", "b"), A = "a", B = "b")
  old <- branch_rate_pipeline(maps, groups, n_points = 3)
  expect_null(old$phenomic)
  integrated <- branch_rate_pipeline(maps, groups, n_points = 3,
                                      phenomic = list(groups = c("A", "B")))
  separate <- phenomic_metrics(old, groups = c("A", "B"), n_points = 3)
  expect_equal(integrated$phenomic, separate)
  integrated$phenomic <- NULL
  expect_equal(integrated, old)
  branch_only <- branch_rate_pipeline(maps, groups, marginal = FALSE,
                                      phenomic = list(groups = "A"))
  expect_null(branch_only$marginal_summary)
  expect_true(is.list(branch_only$phenomic$samples))
  expect_error(branch_rate_pipeline(maps, phenomic = TRUE), "named options")
  expect_error(branch_rate_pipeline(maps, phenomic = list(time_points = 0)), "named options")
  expect_error(branch_rate_pipeline(maps, phenomic = list(groups = "A", groups = "A")),
                "named options")
})

test_that("fitted maps and disk counts give identical phenomic results", {
  tree <- asr_test_tree()
  fitted <- fit_models(tree, asr_test_data(), models = "ER", root_p = "equal",
                       nstarts = 0, seed = 901, on_error = "stop")
  mapped <- paramo_stm(fitted, n_stm = 4, batch_size = 2, store = "both",
                       seed = 902, on_error = "stop")
  groups <- list(A = "first", B = "second")
  raw <- branch_rate_pipeline(read_stm(mapped), groups, n_points = 4, phenomic = list())
  direct <- branch_rate_pipeline(mapped, groups, n_points = 4, phenomic = list())
  expect_equal(raw$phenomic, direct$phenomic)
  directory <- tempfile("phenomic-counts-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  disk <- paramo_stm(fitted, n_stm = 4, batch_size = 2, store = "branch_counts",
                     outdir = directory, seed = 902, on_error = "stop")
  counted <- branch_rate_pipeline(disk, groups, n_points = 4, phenomic = list())
  expect_equal(raw$phenomic, counted$phenomic)
  normalized <- branch_rate_pipeline(mapped, groups, normalize = "characters",
                                      n_points = 4, phenomic = list())
  expect_equal(raw$phenomic$samples, normalized$phenomic$samples)
})
