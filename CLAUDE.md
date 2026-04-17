# CLAUDE.md — TMESplit R package

## Git commit rules

**Hard constraint: do NOT add `Co-Authored-By: Claude ... <noreply@anthropic.com>` (or any Anthropic/Claude co-author trailer) to commit messages.** This applies to direct commits and to commits made by any dispatched subagent. When dispatching subagents that will commit, explicitly include this rule in their prompt — subagents may otherwise inherit the default Claude Code commit template that includes the trailer.

If a commit accidentally contains the trailer and has not been pushed, rewrite history (`git rebase -i` or `git filter-branch`) to strip it before pushing.

## Package layout

- S4 Bioconductor-style R package wrapping the Python `tmesplitpy` package via `basilisk`/`reticulate`.
- Main class: `TMESplitResult` (slots: W_shared, W_specific, H_fractions, programs, p_value, p_value_per_group, groups, k_shared, k_specific, call, metadata).
- Main entrypoint: `tmesplit()` with methods for SummarizedExperiment / SingleCellExperiment / matrix / data.frame — all converge to internal `.run_tmesplit()` which calls Python then hydrates via `.hydrate_result()`.
- Plots: `plotPrograms` (ComplexHeatmap), `plotActivities`/`plotSignificance`/`plotNetwork` (ggplot2). All return `invisible()` so they're pipeline-composable but also print to the active device for interactive RStudio use.
- CI runs mocked unit tests only (no basilisk); full basilisk roundtrip tests live in `tests/testthat/test-integration.R` and are guarded by `skip_on_ci()`.
