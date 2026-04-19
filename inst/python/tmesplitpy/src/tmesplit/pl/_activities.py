"""Per-patient program activity (H_fractions) by group."""

from __future__ import annotations

import numpy as np
import pandas as pd

from tmesplit.pl._common import _require_fitted, _require_matplotlib


def activities(model, figsize: tuple | None = None, ax=None):
    """Boxplot of H_fractions (program activity per patient), split by group."""
    _require_matplotlib()
    _require_fitted(model)

    import matplotlib.pyplot as plt
    import seaborn as sns

    H = model.H_fractions
    groups = model.groups
    group_key = model._group_key
    labels = model.adata.obs[group_key].astype(str).to_numpy()

    rows = []
    for g in groups:
        mask = labels == g
        Hg = H[g] if isinstance(H, dict) else np.asarray(H)[mask]
        for patient_idx in range(Hg.shape[0]):
            for prog_idx in range(Hg.shape[1]):
                rows.append({
                    "group": g,
                    "program": f"p{prog_idx}",
                    "activity": Hg[patient_idx, prog_idx],
                })
    df = pd.DataFrame(rows)

    if figsize is None:
        n_progs = df["program"].nunique()
        figsize = (max(6, 0.8 * n_progs + 2), 4.5)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    sns.boxplot(data=df, x="program", y="activity", hue="group", ax=ax)
    ax.set_ylabel("program activity")
    ax.set_title("Program activity per patient")
    return ax
