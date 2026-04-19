# TMESplit 0.99.1

* Vendor the `tmesplitpy` Python source into `inst/python/tmesplitpy/`
  and install it via `basilisk::BasiliskEnvironment(paths = ...)`
  instead of a `pip` git URL. The old `pip = "tmesplitpy @ git+..."`
  spec was rejected by `basilisk:::.check_versions()` (which requires
  every `pip` entry be pinned with `==version`), so no end user could
  provision the env. Will switch back to a pinned PyPI release once
  `tmesplitpy` is published there.

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
