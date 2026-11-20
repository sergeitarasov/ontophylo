library(ontophylo)
library(ggplot2)
source('../branch-rates/figure-style.R')
legacy <- new.env(parent = baseenv())
legacy$quantile <- stats::quantile
source('legacy-reference.R', local = legacy)
set.seed(88120)
for (i in 1:100) {
  rates <- matrix(rexp(120), 30, 4, dimnames = list(NULL, LETTERS[1:4]))
  for (direction in c('over', 'under')) {
    old <- legacy$get_enrichment_mask(rates, 0.05, 0.95, direction)
    current <- get_enrichment_mask(rates, 0.05, 0.95, direction)
    stopifnot(all(old == current))
  }
  patterns <- matrix(rbinom(120, 1, 0.3), 30, 4)
  states <- character(30)
  for (b in 1:30) states[b] <- paste0(patterns[b, ], collapse = '')
  stopifnot(abs(legacy$diversity_metric(states, 'Neff_scaled') -
                  phenomic_diversity(patterns)$Neff) < 1e-12)
  stopifnot(abs(legacy$relative_null_metric(states, 'scaled') -
                  phenomic_coordination(patterns)$Pc) < 1e-12)
}
writeLines('100 fixed-seed cases passed: both directional masks, Neff and scaled Se/Pc match the original formulas.',
             'results/legacy-compatibility.txt')
scenarios <- c('homogeneous', 'coordinated', 'mosaic', 'decelerated')
metric_names <- c('Neff_plus', 'Neff_minus', 'Pc')
metric_labels <- c('Neff+', 'Neff-', 'Pc')
labels <- c('Homogeneous', 'Coordinated', 'Mosaic', 'Decelerated')
all_tables <- list()
index <- 0L
panels <- character()
sensitivity <- list()
for (i in seq_along(scenarios)) {
  x <- readRDS(paste0('results/', scenarios[i], '.rds'))
  # Recalculate with the final implementation; cached repeated slices must agree.
  updated <- phenomic_metrics(x$inferred, n_points = 61)
  stopifnot(isTRUE(all.equal(updated$samples, x$inferred$phenomic$samples)))
  strict <- phenomic_metrics(x$inferred, n_points = 61, ties = 'exclude')
  strict_low <- mean(strict$regimes$decelerated)
  inclusive_low <- mean(updated$regimes$decelerated)
  sensitivity[[i]] <- data.frame(scenario = scenarios[i],
                                 inclusive_slow_fraction = inclusive_low,
                                 strict_slow_fraction = strict_low)
  for (j in seq_along(metric_names)) {
    panel <- paste(labels[i], metric_labels[j], sep = ': ')
    panels <- c(panels, panel)
    for (stage in c('generating', 'realized', 'inferred')) {
      result <- x[[stage]]
      if (stage != 'generating') result <- result$phenomic
      tab <- result$summary[result$summary$metric == metric_names[j], ]
      index <- index + 1L
      all_tables[[index]] <- data.frame(scenario = scenarios[i], panel = panel,
                                         stage = stage, tab, row.names = NULL)
    }
  }
}
curves <- do.call(rbind, all_tables)
curves$panel <- factor(curves$panel, levels = panels)
curves$stage <- factor(curves$stage, levels = c('generating', 'realized', 'inferred'))
write.csv(curves, 'results/metric-curves.csv', row.names = FALSE)
write.csv(do.call(rbind, sensitivity), 'results/tie-sensitivity.csv', row.names = FALSE)
write.csv(curves[curves$time == 0, ], 'results/present-metrics.csv', row.names = FALSE)
# Duplicate coordinates make the interval band a step function too.
bands <- curves[curves$stage == 'inferred', ]
step_bands <- list()
for (i in seq_along(panels)) {
  tab <- bands[bands$panel == panels[i], ]
  tab <- tab[order(tab$time, decreasing = TRUE), ]
  n <- nrow(tab)
  old <- tab[c(rep(seq_len(n - 1), each = 2), n), ]
  old$time <- c(tab$time[1], rep(tab$time[-1], each = 2))
  step_bands[[i]] <- old
}
bands <- do.call(rbind, step_bands)
figure <- ggplot(curves, aes(time, mean, colour = stage, linetype = stage)) +
  geom_ribbon(data = bands, aes(x = time, ymin = lower, ymax = upper),
                inherit.aes = FALSE, fill = '#0072B2', alpha = 0.17) +
  geom_step(linewidth = 0.4) +
  facet_wrap(~panel, ncol = 3, scales = 'free_y') +
  scale_x_reverse(breaks = c(150, 75, 0), expand = expansion(mult = 0.02)) +
  scale_colour_manual(values = c(generating = 'black', realized = '#E69F00', inferred = '#0072B2'),
                        labels = c('Generating rates', 'True histories', 'Inferred histories')) +
  scale_linetype_manual(values = c(generating = 'dotted', realized = 'dashed', inferred = 'solid'),
                          labels = c('Generating rates', 'True histories', 'Inferred histories')) +
  labs(x = 'Age (Ma)', y = 'Metric value', colour = NULL, linetype = NULL) +
  theme_nature() + theme(legend.position = 'bottom', strip.text = element_text(size = 9),
                          panel.spacing = grid::unit(3, 'mm'))
save_nature(figure, 'results/phenomic-pilot.pdf', width_mm = 183, height_mm = 170)
caption <- paste('Phenomic-metric smoke tests on one fixed 100-tip tree.',
  'Each of four prespecified datasets contains two groups of ten binary characters.',
  'ER models were fitted separately with equal root frequencies, followed by 25 maps per character.',
  'The 5th/95th-percentile rule includes ties and leaves constant columns unclassified.',
  'Black dotted curves apply that same rule to the generating rate vectors; orange dashed curves use the actual simulated histories;',
  'blue curves and bands summarize inferred histories and their pointwise 95% conditional mapping quantiles.',
  'Generating-rate patterns, realized-count patterns and inferred-count patterns are different targets.',
  'Their disagreement includes sampling variation in event counts and the quantile rule, not just reconstruction error.',
  'The four single datasets are smoke tests, not power or coverage estimates. Axes are linear.')
writeLines(caption, 'results/figure-caption.txt')
