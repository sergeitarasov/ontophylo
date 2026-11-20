test_that("ER likelihood and MLE agree with independent exact enumeration", {
  tree <- asr_test_tree()
  x <- asr_test_data(tree)[, 1, drop = FALSE]
  fit <- fit_models(tree, x, models = c("ER", "SYM", "ARD"),
                    root_p = "equal", nstarts = 0, ancestral = "marginal")
  one <- fit$characters$first
  objective <- function(rate) {
    Q <- matrix(c(-rate, rate, rate, -rate), 2, 2)
    asr_exact(tree, one$tips_prior, Q, "equal")$logLik
  }
  direct <- optimize(objective, c(1e-5, 2), maximum = TRUE, tol = 1e-8)
  expect_equal(one$Q_matrix[1, 2], direct$maximum, tolerance = 1e-4)
  expect_equal(one$model_table$logLik[1], direct$objective, tolerance = 1e-7)
  expect_equal(one$model_table$AICc[1], -2 * direct$objective + 2 + 4/6,
               tolerance = 1e-7)
  exact <- asr_exact(tree, one$tips_prior, one$Q_matrix, "equal")
  expect_equal(unname(one$ancestral), exact$marginal, tolerance = 1e-6)
  expect_identical(one$selected_model, "ER")
  expect_equal(one$model_table$equivalent_to[2], "ER")
  expect_equal(sum(one$model_table$weight, na.rm = TRUE), 1)
  expect_true(one$comparison_complete)
})

test_that("ARD root treatments and ambiguous states agree with exact likelihoods", {
  tree <- ape::read.tree(text = "((a:0.8,b:0.7):0.4,(c:0.5,d:0.9):0.6);")
  data <- data.frame(trait = c("blue", "blue&red", "red", "?"),
                      row.names = tree$tip.label)
  roots <- list("equal", "yang", "maddfitz", c(blue = 0.8, red = 0.2))
  for (root in roots) {
    fits <- fit_models(tree, data, models = "ARD", root_p = root,
                       criterion = "AIC", nstarts = 0, ancestral = "marginal",
                       control = list(upper.bound = 2))
    one <- fits$characters$trait
    exact <- asr_exact(tree, one$tips_prior, one$Q_matrix, root)
    expect_equal(one$model_table$logLik, exact$logLik, tolerance = 1e-6)
    expect_equal(unname(one$ancestral), exact$marginal, tolerance = 1e-5)
    expect_equal(rowSums(one$Q_matrix), c(blue = 0, red = 0), tolerance = 1e-12)
    expect_equal(one$tips_prior["b", ], c(blue = 1, red = 1))
  }
})

test_that("named custom models preserve forbidden transitions and rate constraints", {
  tree <- asr_test_tree()
  x <- data.frame(trait = c("red", "red", "green", "green", "blue", "blue", "green", "red"),
                   row.names = tree$tip.label)
  index <- matrix(c(0, 1, 0, 1, 0, 2, 0, 2, 0), 3, 3,
                    dimnames = list(c("red", "green", "blue"), c("red", "green", "blue")))
  fits <- fit_models(tree, x, models = "ordered", custom_models = list(ordered = index),
                     root_p = "equal", criterion = "AIC", nstarts = 0)
  one <- fits$characters$trait
  Q <- one$Q_matrix
  expect_equal(Q["red", "blue"], 0)
  expect_equal(Q["blue", "red"], 0)
  expect_equal(Q["red", "green"], Q["green", "red"])
  expect_equal(Q["blue", "green"], Q["green", "blue"])
  expect_equal(one$model_table$npar, 2L)
  # Three states on seven internal nodes: 3^7 independent configurations.
  exact <- asr_exact(tree, one$tips_prior, Q, "equal")
  expect_equal(one$model_table$logLik, exact$logLik, tolerance = 1e-6)
})

test_that("data order, character overrides and root names are respected", {
  tree <- asr_test_tree()
  data <- asr_test_data(tree)
  a <- fit_models(tree, data, models = "ER", nstarts = 1, seed = 33)
  b <- fit_models(tree, data[8:1, 2:1], models = "ER", nstarts = 1, seed = 33)
  expect_equal(a$characters$first$Q_matrix, b$characters$first$Q_matrix)
  expect_equal(a$characters$second$Q_matrix, b$characters$second$Q_matrix)
  c <- fit_models(tree, data, models = "ER", nstarts = 0,
                  model_overrides = list(second = "ARD"),
                  root_p = list(first = "equal", second = c("1" = .8, "0" = .2)))
  expect_equal(c$characters$second$root_p, c("0" = .2, "1" = .8))
  expect_equal(c$characters$second$selected_model, "ARD")
  expect_error(fit_models(tree, data[-1, ]), "match tree")
  expect_error(fit_models(tree, data, control = list(typo = TRUE)), "control")
  expect_error(fit_models(tree, data, models = "bogus"), "Invalid")
})

test_that("failures, constant characters and undefined AICc remain visible", {
  tree <- asr_test_tree()
  data <- asr_test_data(tree)
  data$constant <- "0"
  data$missing <- NA_character_
  data$inapplicable <- "-"
  expect_warning(fits <- fit_models(tree, data, models = "ER", nstarts = 0),
                  "failed or were skipped")
  expect_equal(fits$status$status, c("ok", "ok", "skipped", "skipped", "failed"))
  expect_match(fits$characters$missing$reason, "missing")
  expect_error(fit_models(tree, data, models = "ER", nstarts = 0, on_error = "stop"),
                 "failed or were skipped")
  small <- ape::read.tree(text = "(a:1,b:1);")
  expect_warning(bad <- fit_models(small, c(a = "0", b = "1"), models = "ER", nstarts = 0),
                  "failed or were skipped")
  expect_equal(bad$status$status, "failed")
  expect_match(bad$status$reason, "No valid model under AICc")
  expect_true(is.na(bad$characters[[1]]$model_table$AICc))
  tree$edge.length[1] <- 1e-6
  expect_error(fit_models(tree, data), "Rescale")
})

test_that("fit checkpoints, parallelism and RNG are reproducible", {
  tree <- asr_test_tree()
  data <- asr_test_data(tree)
  out <- tempfile("fits-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  set.seed(991)
  old <- .Random.seed
  serial <- fit_models(tree, data, models = c("ER", "ARD"), nstarts = 1,
                       seed = 71, outdir = out)
  expect_identical(.Random.seed, old)
  resumed <- fit_models(tree, data, models = c("ER", "ARD"), nstarts = 1,
                         seed = 71, outdir = out, resume = TRUE)
  expect_identical(serial, resumed)
  expect_error(fit_models(tree, data, models = "ER", outdir = out, resume = TRUE),
                 "checkpoint mismatch")
  parallel <- fit_models(tree, data, models = c("ER", "ARD"), nstarts = 1,
                          seed = 71, workers = 2)
  expect_equal(serial$characters, parallel$characters)
  expect_identical(.Random.seed, old)
})

test_that("a failed candidate cannot silently make a complete comparison", {
  tree <- asr_test_tree()
  wrong <- matrix(1, 3, 3)
  expect_warning(fits <- fit_models(tree, asr_test_data(tree)[, 1, drop = FALSE],
                    models = c("ER", "wrong"), custom_models = list(wrong = wrong),
                    root_p = "equal", nstarts = 0), "Incomplete model comparisons")
  expect_false(fits$characters$first$comparison_complete)
  expect_equal(fits$characters$first$model_table$status, c("ok", "failed"))
  expect_match(fits$characters$first$model_table$error[2], "state names")
})
