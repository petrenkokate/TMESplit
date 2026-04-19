"""Cell-type loading heatmap for shared and group-specific programs."""

from __future__ import annotations

import numpy as np

from tmesplit.pl._common import _require_fitted, _require_matplotlib


def programs(model, figsize: tuple | None = None, cmap: str = "magma",
             ax=None):
    """Heatmap of W matrices: rows = cell types, columns = programs.

    Columns are annotated with "shared" or a group label.
    """
    _require_matplotlib()
    _require_fitted(model)

    import matplotlib.pyplot as plt
    import seaborn as sns

    W_s = model.W_shared
    W_spec = model.W_specific
    cell_types = list(model.adata.var.index)

    columns = [W_s] if W_s.size else []
    col_labels = [f"shared_{i}" for i in range(W_s.shape[1])] if W_s.size else []
    for g, Wg in W_spec.items():
        if Wg.size == 0:
            continue
        columns.append(Wg)
        col_labels.extend([f"{g}_{i}" for i in range(Wg.shape[1])])

    if not columns:
        raise ValueError("No programs to plot (neither shared nor specific).")
    W = np.concatenate(columns, axis=1)

    if figsize is None:
        figsize = (max(6, 0.5 * W.shape[1] + 3), max(6, 0.25 * W.shape[0]))
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    sns.heatmap(W, yticklabels=cell_types, xticklabels=col_labels,
                cmap=cmap, ax=ax, cbar_kws={"label": "loading"})
    ax.set_xlabel("program")
    ax.set_ylabel("cell type")
    ax.set_title("TMESplit programs")
    return ax
