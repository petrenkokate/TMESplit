# TMESplit

Formal testing of group-specific multicellular coordination programs from
cell-type composition data. Bioconductor-style R interface; wraps the
[`TMESplitpy`](https://github.com/petrenkokate/TMESplitpy) Python
implementation via [basilisk](https://bioconductor.org/packages/basilisk).

**Status:** pre-release (0.99.0). Scaffold only — see `task_plan.md` in the
EcoSplit paper repo for roadmap.

## Installation

```r
# not yet released
# remotes::install_github("petrenkokate/TMESplit")
```

On first use, `basilisk` provisions a sandboxed conda env containing
`tmesplitpy`. No user-side Python setup required.

## Minimal usage (planned API)

```r
library(TMESplit)
library(SummarizedExperiment)

res <- tmesplit(se, group_col = "condition",
                beta_loss = 0.5, n_perm = 2000)

res                     # S4 summary
programs(res)           # DataFrame of programs + FDR
pvalue(res, level = "global")
plotPrograms(res)
```

## License

MIT. See [LICENSE.md](LICENSE.md).
