# TMESplitpy

**T**umor **M**icro**E**nvironment **Split** — formal testing of group-specific multicellular coordination programs from cell-type compositions.

TMESplit detects and tests whether the coordination structure of cell-type compositions differs between biological groups (e.g. tumor vs normal, responders vs non-responders, cancer subtypes). It decomposes compositions into *shared* and *group-specific* programs via per-group NMF, cross-correlation matching, and a permutation test with three-level hierarchical significance reporting.

This is the Python implementation. An R/Bioconductor companion lives at [TMESplit](https://github.com/petrenkokate/TMESplit).

## Status

Alpha. API is not yet frozen. Paper in preparation.

## Installation

```bash
pip install tmesplitpy            # once published to PyPI
```

From source:

```bash
git clone git@github.com:petrenkokate/TMESplitpy.git
cd TMESplitpy
pip install -e ".[dev,plot]"
```

## Quickstart

```python
import scanpy as sc
import tmesplit as tme

adata = sc.read_h5ad("cells.h5ad")

comp = tme.pp.aggregate_to_composition(
    adata, patient_key="patient_id", celltype_key="celltype_L2"
)

tme.TMESplit.setup_anndata(comp, group_key="condition")
model = tme.TMESplit(comp, n_perm=2000, random_state=42)
model.train(n_jobs=8)

model.programs        # tidy DataFrame of programs + FDR
model.p_value         # global (Level 1) p-value
tme.pl.programs(model)
```

## License

MIT — see [LICENSE](LICENSE).

## Citation

Paper in preparation. Will update on preprint release.
