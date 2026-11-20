test_that("mapping preserves observations, branches and conditional rate samples", {
  tree <- asr_test_tree()
  fits <- fit_models(tree, asr_test_data(tree), models = "ER", root_p = "equal", nstarts = 0)
  mapped <- paramo_stm(fits, n_stm = 8, batch_size = 3, seed = 51, store = "both")
  maps <- read_stm(mapped)
  counts <- read_stm(mapped, type = "branch_counts")
  for (id in names(maps)) {
    for (i in seq_along(maps[[id]])) {
      map <- maps[[id]][[i]]
      expect_identical(map$edge, fits$tree$edge)
      expect_equal(counts[[id]][, i], count_branch_changes(map), ignore_attr = TRUE)
      for (b in seq_len(nrow(tree$edge))) {
        expect_equal(sum(map$maps[[b]]), map$edge.length[b])
        tip <- map$edge[b, 2]
        if (tip <= length(tree$tip.label)) {
          state <- tail(names(map$maps[[b]]), 1)
          expect_equal(fits$characters[[id]]$tips_prior[tip, state], 1)
        }
      }
    }
  }
  combined <- branch_rates(maps)
  expected <- (counts$first + counts$second) / fits$tree$edge.length
  expect_equal(combined$rate_samples, expected)
  expect_equal(branch_rates(mapped), combined)
})

test_that("conditional root frequencies agree with exact ancestral probabilities", {
  tree <- ape::read.tree(text = "(a:0.8,b:1.2);")
  Q <- matrix(c(-.3, .8, .3, -.8), 2, 2,
                dimnames = list(c("0", "1"), c("0", "1")))
  tip <- matrix(c(1, 0, 0, 1), 2, 2, dimnames = list(c("a", "b"), c("0", "1")))
  roots <- list("equal", "yang", "maddfitz", c("0" = .8, "1" = .2))
  for (root in roots) {
    fits <- asr_known_fit(tree, Q, tip, root)
    mapped <- paramo_stm(fits, n_stm = 1200, batch_size = 400, seed = 82)
    maps <- read_stm(mapped)$known
    root0 <- logical(length(maps))
    for (i in seq_along(maps)) root0[i] <- names(maps[[i]]$maps[[1]])[1] == "0"
    expected <- asr_exact(tree, tip, Q, root)$marginal[1, 1]
    se <- sqrt(expected * (1 - expected) / length(maps))
    expect_lt(abs(mean(root0) - expected), 6 * se + 0.005)
  }
})

test_that("mapped change counts obey the analytical binary bridge distribution", {
  tree <- ape::read.tree(text = "(a:1,b:1);")
  rate <- .7
  Q <- matrix(c(-rate, rate, rate, -rate), 2, 2,
                dimnames = list(c("0", "1"), c("0", "1")))
  tip <- matrix(c(1, 0, 0, 1), 2, 2, dimnames = list(c("a", "b"), c("0", "1")))
  fits <- asr_known_fit(tree, Q, tip, c("0" = 1, "1" = 0))
  mapped <- paramo_stm(fits, n_stm = 2000, batch_size = 500, seed = 211,
                        store = "branch_counts")
  counts <- read_stm(mapped, type = "branch_counts")$known
  # Symmetric binary CTMC: Poisson(q*t), conditioned on even/odd jump counts.
  expected <- c(rate * tanh(rate), rate / tanh(rate))
  variance <- rate^2 + expected - expected^2
  expect_true(all(counts[1, ] %% 2 == 0))
  expect_true(all(counts[2, ] %% 2 == 1))
  expect_true(all(abs(rowMeans(counts) - expected) < 6 * sqrt(variance/2000)))
  expect_error(read_stm(mapped), "not stored")
  rates <- branch_rates(mapped, marginal = FALSE)
  expect_equal(unname(rates$rate_samples), unname(counts))
})

test_that("batch sizes, parallel workers and resumption retain identical histories", {
  tree <- asr_test_tree()
  fits <- fit_models(tree, asr_test_data(tree), models = "ER", nstarts = 0)
  set.seed(121)
  old <- .Random.seed
  a <- paramo_stm(fits, n_stm = 9, batch_size = 9, seed = 17, store = "both")
  b <- paramo_stm(fits, n_stm = 9, batch_size = 2, seed = 17, store = "both", workers = 2)
  expect_identical(read_stm(a), read_stm(b))
  expect_identical(.Random.seed, old)
  out <- tempfile("maps-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  partial <- paramo_stm(fits, n_stm = 4, batch_size = 2, seed = 17,
                         store = "both", outdir = out)
  resumed <- paramo_stm(fits, n_stm = 9, batch_size = 3, seed = 17,
                         store = "both", outdir = out, resume = TRUE)
  expect_identical(read_stm(a), read_stm(resumed))
  expect_identical(read_stm(a, type = "branch_counts"), read_stm(resumed, type = "branch_counts"))
  expect_true(is.null(resumed$characters$first$maps))
  expect_equal(branch_rates(resumed), branch_rates(a))
  expect_warning(changed <- paramo_stm(fits, n_stm = 9, seed = 18,
                      store = "both", outdir = out, resume = TRUE), "failed")
  expect_true(all(changed$status$status == "failed"))
  expect_true(all(grepl("mismatch", changed$status$reason)))
  # Simulate interruption by removing the final complete chunk of one character.
  unlink(tail(resumed$characters$first$files, 1))
  recovered <- paramo_stm(fits, n_stm = 9, batch_size = 1, seed = 17,
                           store = "both", outdir = out, resume = TRUE)
  expect_identical(read_stm(a), read_stm(recovered))
})

test_that("forbidden transitions remain absent from sampled histories", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  Q <- matrix(c(-.5, .5, 0, .5, -1, .5, 0, .5, -.5), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  tip <- rbind(a = c(1, 0, 0), b = c(0, 1, 0), c = c(0, 0, 1), d = c(0, 0, 1))
  colnames(tip) <- colnames(Q)
  fits <- asr_known_fit(tree, Q, tip)
  maps <- read_stm(paramo_stm(fits, n_stm = 40))$known
  for (map in maps) {
    for (branch in map$maps) {
      states <- names(branch)
      if (length(states) > 1) {
        transitions <- cbind(states[-length(states)], states[-1])
        expect_true(all(Q[transitions] > 0))
      }
    }
  }
})
