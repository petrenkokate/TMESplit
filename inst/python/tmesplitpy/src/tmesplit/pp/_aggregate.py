"""Aggregate cell-level AnnData to patient x cell-type compositions."""

from __future__ import annotations

from collections.abc import Sequence

import anndata as ad
import numpy as np
import pandas as pd


def aggregate_to_composition(
    adata: ad.AnnData,
    patient_key: str,
    celltype_key: str,
    patient_meta_cols: Sequence[str] | None = None,
) -> ad.AnnData:
    """Aggregate a cell-level AnnData to patient x cell-type compositions.

    Parameters
    ----------
    adata
        Cell-level AnnData. Rows are cells; must have ``patient_key`` and
        ``celltype_key`` as columns in ``adata.obs``.
    patient_key
        Column in ``adata.obs`` that identifies each patient/sample.
    celltype_key
        Column in ``adata.obs`` that identifies each cell's cell type.
    patient_meta_cols
        Columns in ``adata.obs`` to carry over as patient-level metadata.
        Each must be constant within patient. If ``None``, every column
        that is constant within patient is carried over automatically.

    Returns
    -------
    AnnData
        Shape ``(n_patients, n_celltypes)``. ``X`` holds proportions that
        sum to 1 per row. ``obs`` carries patient metadata; ``var`` is
        indexed by cell-type names.
    """
    for k in (patient_key, celltype_key):
        if k not in adata.obs.columns:
            raise KeyError(f"{k!r} not found in adata.obs")

    obs = adata.obs.copy()
    obs[patient_key] = obs[patient_key].astype(str)
    obs[celltype_key] = obs[celltype_key].astype(str)

    counts = (
        obs.groupby([patient_key, celltype_key], observed=True)
        .size()
        .unstack(fill_value=0)
        .astype(float)
    )
    row_sums = counts.sum(axis=1)
    row_sums[row_sums == 0] = 1.0
    props = counts.div(row_sums, axis=0)

    if patient_meta_cols is None:
        candidates = [c for c in obs.columns if c != celltype_key]
        nunique = obs.groupby(patient_key, observed=True)[candidates].nunique(dropna=False)
        keep = [c for c in candidates if (nunique[c] <= 1).all()]
    else:
        keep = [c for c in patient_meta_cols if c in obs.columns and c not in (celltype_key,)]
        if patient_key not in keep:
            keep = [patient_key, *keep]

    patient_meta = (
        obs[keep].drop_duplicates(subset=[patient_key]).set_index(patient_key)
    )
    patient_meta = patient_meta.loc[props.index]
    if patient_key in patient_meta.columns:
        patient_meta = patient_meta.drop(columns=[patient_key])

    var = pd.DataFrame(index=props.columns)
    var.index.name = celltype_key

    return ad.AnnData(
        X=props.values.astype(np.float64),
        obs=patient_meta,
        var=var,
    )
