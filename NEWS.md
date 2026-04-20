# TMESplit 0.99.2

* Plotting refresh (Phase 3.7):
    * New `tme_palette` / `tme_diverging` color constants re-used across all
      plot functions.
    * `plotPrograms()` now supports `scale`, `top_n`, `transpose`,
      `program_names`, and `cluster_rows` / `cluster_columns`; the
      diverging colour scale matches the `sex_diff/6_TIME.Rmd` reference.
    * `plotActivities()` annotates each program with a bracket + p-value.
      Shared programs use a between-group test on H activities
      (`test_method = "wilcoxon"` by default, `"t"` and `"kruskal"` opt-in;
      forced to Kruskal-Wallis when G >= 3). Group-specific programs fall
      back to the Level 3 `combined_p` from `programs_table` with a `"(L3)"`
      marker so the test source is explicit.
    * `plotSignificance()` recoloured against the shared palette.
    * New `plotSampleHeatmap()` — samples x cell-types heatmap split by
      dominant program with user-supplied annotations. Defaults now sort
      cell types by their dominant program (argmax of normalised
      `[W_shared | W_specific]`) with within-program loading descending,
      and skip row/column clustering, producing a diagonal layout; set
      `cluster_samples = TRUE` or `cluster_cell_types = TRUE` to opt back
      into Ward clustering. New `row_split = TRUE` default slices the rows
      by dominant program (mirrors the existing column split).
    * `plotNetwork()` replaced with an `igraph` Fruchterman-Reingold
      force-directed layout; nodes are coloured by dominant program so
      cell types from the same program cluster spatially. Default
      `threshold` lowered from `0.25` to `0.1` — real-data TME programs are
      often diffuse, and the old default excluded most secondary
      contributors. The docstring now spells out the co-loading edge
      semantics.
* `TMESplitResult@H_fractions` is now a named list of per-group matrices
  (was a padded combined matrix, which misrepresented specific-program
  column labels). The raw composition is stored under
  `metadata(result)$composition` for downstream plots.
* `DESCRIPTION` gains `igraph` as an import (force-directed layout).
* Vignette `TMESplit-analysis.Rmd` refreshed with the new API
  (`plotPrograms(top_n=)`, `plotActivities(test_method=)`, a new
  `plotSampleHeatmap` section with annotations, and the co-loading
  description for `plotNetwork`). Mock test fixture switched to the new
  per-group `H_fractions` list layout and populated with
  `metadata$composition` so `plotSampleHeatmap` is covered in unit tests.
* Bug fix: `plotPrograms(transpose=TRUE)` no longer errors with "formal
  argument 'which' matched by multiple actual arguments" (duplicate
  `which = "row"` passed to `ComplexHeatmap::rowAnnotation`).

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
