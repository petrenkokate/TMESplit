"""Cell-type co-activation network derived from program loadings."""

from __future__ import annotations

import numpy as np

from tmesplit.pl._common import _require_fitted, _require_matplotlib


def network(model, top_k: int = 3, threshold: float = 0.3,
            figsize: tuple | None = None, ax=None):
    """Cell-type co-loading graph: edge if two cell types co-load (> threshold)
    on the same top-``top_k`` component of any program.

    Minimal matplotlib-only implementation (no networkx dep).
    """
    _require_matplotlib()
    _require_fitted(model)

    import matplotlib.pyplot as plt

    cell_types = list(model.adata.var.index)
    W_parts = [model.W_shared] if model.W_shared.size else []
    for Wg in model.W_specific.values():
        if Wg.size:
            W_parts.append(Wg)
    if not W_parts:
        raise ValueError("No programs to plot.")
    W = np.concatenate(W_parts, axis=1)

    Wn = W / (np.linalg.norm(W, axis=0, keepdims=True) + 1e-12)
    edges = {}
    for j in range(Wn.shape[1]):
        top = np.argsort(Wn[:, j])[-top_k:]
        for a in top:
            for b in top:
                if a >= b or Wn[a, j] < threshold or Wn[b, j] < threshold:
                    continue
                key = (int(a), int(b))
                edges[key] = max(edges.get(key, 0.0), float(min(Wn[a, j], Wn[b, j])))

    n = len(cell_types)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    xs, ys = np.cos(angles), np.sin(angles)

    if figsize is None:
        figsize = (7, 7)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    for (a, b), w in edges.items():
        ax.plot([xs[a], xs[b]], [ys[a], ys[b]],
                color="gray", alpha=0.3 + 0.7 * w, linewidth=0.5 + 3 * w)
    ax.scatter(xs, ys, s=80, color="#e74c3c", zorder=3)
    for i, name in enumerate(cell_types):
        ax.text(xs[i] * 1.08, ys[i] * 1.08, name, ha="center", va="center", fontsize=8)

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Cell-type co-loading network")
    return ax
