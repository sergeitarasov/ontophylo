# Hand-counted histories give exact expected branch rates.
make_test_map <- function(counts, tree = ape::read.tree(text = "(a:2,b:2);")) {
  tree$maps <- vector("list", nrow(tree$edge))
  for (b in seq_along(counts)) {
    n_segments <- counts[b] + 1L
    states <- rep(c("0", "1"), length.out = n_segments)
    durations <- rep(tree$edge.length[b] / n_segments, n_segments)
    names(durations) <- states
    tree$maps[[b]] <- durations
  }
  class(tree) <- c("simmap", "phylo")
  tree
}

test_that("original character counts are summed within each replicate", {
  chars <- list(
    first = list(make_test_map(c(0, 2)), make_test_map(c(2, 0))),
    second = list(make_test_map(c(2, 0)), make_test_map(c(0, 4))))
  expected <- list(c(1, 1), c(1, 2))
  expect_equal(ama_branches_multi_char(chars), expected)
  expect_equal(get_branch_rate_sample_smaps(chars$first), list(c(0, 1), c(1, 0)))
  expect_equal(summarize_branch_rates_multi(chars), c(1, 1.5))
  expect_equal(summarize_branch_rates(chars$first), c(0.5, 0.5))
  expect_equal(summarize_branch_rates_multi(chars, "median"), c(1, 1.5))
  expect_equal(summarize_branch_rates_multi(chars, "sd"), c(0, sqrt(0.5)))
  expect_equal(unname(summarize_branch_rates_multi(chars, "quantile")),
               rbind(c(1, 1.025), c(1, 1.5), c(1, 1.975)))
  expect_error(summarize_branch_rates_multi(chars, "quantile", probs = 2), "probs")

  result <- branch_rates(chars, time_points = c(2, 1, 0), time_scale = "age")
  expect_equal(unname(result$rate_samples), cbind(c(1, 1), c(1, 2)))
  expect_equal(result$branch_summary$mean, c(1, 1.5))
  expect_equal(unname(result$marginal_samples), cbind(rep(1, 3), rep(1.5, 3)))
  expect_equal(result$marginal_summary$mean, rep(1.25, 3))
  expect_equal(result$marginal_summary$n_branches, rep(2L, 3))
  expect_equal(result$characters, names(chars))
  expect_equal(result$n_maps, 2L)
})

test_that("state changes exclude repeated labels and retain reversals", {
  map <- make_test_map(c(2, 0))
  expect_equal(count_branch_changes(map), c(2, 0))
  map$maps[[1]] <- c("0" = 0.5, "0" = 0.5, "1" = 1)
  expect_equal(count_branch_changes(map), c(1, 0))
})

test_that("one-character and one-map cases retain their dimensions", {
  map <- make_test_map(c(2, 0))
  expect_equal(ama_branches_multi_char(list(one = map)), list(c(1, 0)))
  result <- branch_rates(list(one = map))
  expect_equal(dim(result$rate_samples), c(2L, 1L))
  expect_true(all(is.na(result$branch_summary$sd)))
  expect_equal(result$branch_summary$lower, result$branch_summary$upper)
  expect_true(all(is.na(result$marginal_summary$sd)))
  expect_equal(dim(summarize_branch_rates(map, "quantile", 0.5)), c(1L, 2L))
})

test_that("custom time grids use correct ages and nodes count descendants once", {
  tree <- ape::read.tree(text = "((a:2,b:2):8,c:10);")
  rates <- c(1, 3, 5, 7)
  curve <- branch_rate_thru_time(tree, rates, time_points = c(0, 2, 8, 9, 10))
  expect_equal(curve$time, c(10, 8, 2, 1, 0))
  expect_equal(curve$rate, c(4, 4, 5, 5, 5))
  expect_equal(curve$n_branches, c(2L, 2L, 3L, 3L, 3L))
  ages <- branch_rate_thru_time(tree, rates, time_points = c(10, 8, 2, 1, 0),
                                time_scale = "age")
  expect_equal(curve, ages)
  one_time <- branch_rate_thru_time(tree, rates, time_points = 9)
  expect_equal(one_time$time, 1)
  expect_equal(one_time$rate, 5)
  expect_error(branch_rate_thru_time(tree, rates, time_points = 11), "time_points")
  expect_error(branch_rate_thru_time(tree, rates, n_points = 1), "n_points")
  expect_error(branch_rate_thru_time(tree, rates[-1]), "one finite")
})

test_that("marginal uncertainty is computed after averaging each replicate", {
  tree <- ape::read.tree(text = "(a:2,b:2);")
  samples <- list(c(0, 4), c(4, 0))
  result <- branch_rate_thru_time_Multi(tree, samples, n_points = 3,
                                        return_samples = TRUE)
  expect_equal(unname(result$samples), matrix(2, 3, 2))
  expect_equal(result$summary$lower, rep(2, 3))
  expect_equal(result$summary$upper, rep(2, 3))
  expect_equal(result$summary$sd, rep(0, 3))
  old_format <- list(cbind(c(0, 2), c(0, 2)), cbind(c(2, 0), c(2, 0)))
  expect_equal(branch_rate_thru_time_Multi(tree, old_format, n_points = 3),
               result$summary)
})

test_that("incompatible trees, samples and malformed histories fail clearly", {
  map <- make_test_map(c(2, 0))
  expect_error(ama_branches_multi_char(list()), "nonempty list")
  expect_error(ama_branches_multi_char(list(a = list(map), b = list(map, map))),
               "same number")
  expect_error(ama_branches_multi_char(list(a = list(map), a = list(map))), "unique")
  other <- map
  other$edge <- other$edge[2:1, ]
  other$maps <- other$maps[2:1]
  expect_error(ama_branches_multi_char(list(a = list(map), b = list(other))),
               "edge order")
  other <- map
  other$tip.label <- rev(other$tip.label)
  expect_error(ama_branches_multi_char(list(a = list(map), b = list(other))),
               "tip labels")
  other <- map
  other$edge.length[1] <- 3
  expect_error(ama_branches_multi_char(list(a = list(map), b = list(other))),
               "branch lengths")
  other <- map
  other$edge.length[1] <- 0
  expect_error(count_branch_changes(other), "positive lengths")
  other <- map
  other$maps[[1]][1] <- -1
  expect_error(count_branch_changes(other), "positive durations")
  other <- map
  other$maps[[1]][1] <- 5
  expect_error(count_branch_changes(other), "sum to branch lengths")
  other <- map
  names(other$maps[[1]]) <- NULL
  expect_error(count_branch_changes(other), "state labels")
  non_ultrametric <- ape::read.tree(text = "(a:2,b:3);")
  other <- make_test_map(c(2, 0), non_ultrametric)
  expect_error(branch_rates(list(a = list(other))), "ultrametric")
  result <- branch_rates(list(a = list(other)), marginal = FALSE)
  expect_equal(result$branch_summary$mean, c(1, 0))
  expect_null(result$marginal_summary)
  expect_null(result$marginal_samples)
})

test_that("phytools map samples run through the complete workflow", {
  set.seed(191119)
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  Q <- matrix(c(-0.5, 0.5, 0.5, -0.5), 2, 2)
  dimnames(Q) <- list(c("0", "1"), c("0", "1"))
  maps1 <- phytools::make.simmap(tree, c(a = "0", b = "1", c = "0"),
                                Q = Q, pi = c(0.5, 0.5), nsim = 3, message = FALSE)
  maps2 <- phytools::make.simmap(tree, c(a = "1", b = "0", c = "1"),
                                Q = Q, pi = c(0.5, 0.5), nsim = 3, message = FALSE)
  chars <- list(first = maps1, second = maps2)
  result <- branch_rates(chars, n_points = 4)
  # Independent change counts from phytools, across all branches and characters.
  for (i in 1:3) {
    count1 <- phytools::countSimmap(maps1[[i]], message = FALSE)$N
    count2 <- phytools::countSimmap(maps2[[i]], message = FALSE)$N
    reconstructed_count <- sum(result$rate_samples[, i] * tree$edge.length)
    expect_equal(reconstructed_count, count1 + count2)
  }
  # Unit conversion rescales rates and ages, leaving inferred counts unchanged.
  scaled_chars <- chars
  for (j in seq_along(scaled_chars)) {
    for (i in seq_along(scaled_chars[[j]])) {
      scaled_chars[[j]][[i]]$edge.length <- chars[[j]][[i]]$edge.length * 10
      for (b in seq_along(scaled_chars[[j]][[i]]$maps)) {
        scaled_chars[[j]][[i]]$maps[[b]] <- chars[[j]][[i]]$maps[[b]] * 10
      }
    }
  }
  scaled <- branch_rates(scaled_chars, n_points = 4)
  expect_equal(scaled$rate_samples, result$rate_samples / 10)
  expect_equal(scaled$marginal_summary$time, result$marginal_summary$time * 10)
  expect_equal(scaled$marginal_samples, result$marginal_samples / 10)
})
