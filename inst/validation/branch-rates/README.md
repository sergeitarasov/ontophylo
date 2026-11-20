# Branch-specific rate validation

This directory contains a reproducible fixed-tree screening experiment for
`branch_rate_pipeline()`. It tests accounting with true histories separately
from reconstructing heterogeneity using only tip observations.

The exact binary simulator was developed in the Hymenoptera revision workspace
as `simulate-binary-branches.R`. The tree and clade assignment come from the
previously specified seed-186001 experiment: 100 tips, root age 150 Ma, and a
45-tip focal clade. The Newick serialization renumbers internal nodes; the
CSV uses those serialized IDs, with branch durations and regime membership
preserved and checked. A focal branch receives the accelerated rate only when
its entire duration lies below 50 Ma. Branches crossing 50 Ma remain baseline.

Run from this directory with the development ontophylo package installed:

```sh
Rscript --vanilla run-validation.R results 30 50
Rscript --vanilla summarize-validation.R results
```

There are 30 datasets under homogeneous rates and 30 under a 10-fold focal
acceleration, each with 20 binary characters and 50 reconstructed maps per
character. Every character is fitted separately with homogeneous ER Q using
`phytools::fitMk()`. This isolates the downstream implementation from the
package's separate corHMM fitting wrapper. End-to-end tests and the vignette
also exercise `fit_models()` and `paramo_stm()`.

For a two-dataset runtime pilot, set `ONTO_LAST_DATASET=1` before the first
command. Resume with the same output directory and settings. Disjoint ranges
can use `ONTO_FIRST_DATASET` and `ONTO_LAST_DATASET`; seeds depend on scenario,
dataset and character, not scheduling. Run the summarizer after all workers
finish. Do not compare a partial result table with the full validation report.

The saved per-dataset RDS files contain fitted rates, branch and marginal
summaries, generating and realized branch rates, and contrast draws. The CSV
and report separate ranking AUC, inferred fold change, interval coverage, and
the homogeneous-control frequency of excluding equal rates. The figure plots
each dataset, with generating ratios marked independently. AUC is a ranking
measure, not a hypothesis test.

Conditional map intervals do not include stochastic variation across newly
simulated datasets, uncertainty in fitted Q, or tree and sampling uncertainty.
They must not be interpreted as nominal 95% confidence intervals for generating
process parameters. Thirty replicates per condition provide a screening result,
not a precise estimate of tail error rates. No character or dataset was selected
for favourable inference and no thresholds were tuned to the results.

The code retains invariant observations in this benchmark. Unidentified
high-rate fits stop with a saved diagnostic dataset. A stopped run is not
silently resumed by discarding its problematic character.
