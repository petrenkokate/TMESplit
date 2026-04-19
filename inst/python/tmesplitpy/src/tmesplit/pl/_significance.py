"""Hierarchical p-value summary plot (Level 1 / 2 / 3)."""

from __future__ import annotations

import numpy as np

from tmesplit.pl._common import _require_fitted, _require_matplotlib


def significance(model, figsize: tuple | None = None, ax=None):
    """Barplot of -log10(p) across the 3 hierarchical levels."""
    _require_matplotlib()
    _require_fitted(model)

    import matplotlib.pyplot as plt

    labels, values, colors = [], [], []

    labels.append("Global")
    values.append(model.p_value)
    colors.append("#2c3e50")

    for g, p in model.p_value_per_group.items():
        labels.append(f"group: {g}")
        values.append(p)
        colors.append("#3498db")

    programs = model.programs
    if len(programs) > 0 and "combined_p" in programs.columns:
        for _, row in programs.iterrows():
            labels.append(f"{row['group']}/p{int(row['program_idx'])}")
            values.append(float(row["combined_p"]))
            colors.append("#e67e22")

    nlp = [-np.log10(max(float(v), 1e-10)) for v in values]

    if figsize is None:
        figsize = (max(6, 0.6 * len(labels) + 2), 4)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    ax.bar(range(len(labels)), nlp, color=colors)
    ax.axhline(-np.log10(0.05), color="red", linestyle="--", linewidth=0.8,
               label=r"$\alpha=0.05$")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.set_title("Hierarchical significance")
    ax.legend(loc="best", fontsize=8)
    return ax
