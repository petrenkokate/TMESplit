# TMESplit 0.99.0

* Initial pre-release (Bioconductor convention).
* `TMESplitResult` S4 class with slots for shared/group-specific factor
  matrices, patient-level activations, per-program statistics, and
  hierarchical p-values.
* `tmesplit()` generic with methods for `SummarizedExperiment`,
  `SingleCellExperiment`, `matrix`, and `data.frame`.
* `aggregateToComposition()` helper to build a patient x cell-type
  proportions `SummarizedExperiment` from a `SingleCellExperiment`.
* Real basilisk bridge to the Python `tmesplitpy` backend
  (`ecosplit()` + `test_significance()`); results are hydrated into a
  `TMESplitResult` in R.
* Four plotting functions:
  `plotPrograms()` (ComplexHeatmap),
  `plotActivities()` (ggplot2 boxplot),
  `plotSignificance()` (ggplot2 barplot),
  `plotNetwork()` (ggplot2 circular layout).
* Mocked unit tests for CI plus optional full-stack integration tests
  that exercise the basilisk round-trip locally.
* Two vignettes: `TMESplit-quickstart` and `TMESplit-analysis`.
