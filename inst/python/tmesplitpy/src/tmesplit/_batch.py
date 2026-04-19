"""Per-column min-max normalization and per-batch correction methods."""

from __future__ import annotations

import numpy as np


def minmax_normalize(X: np.ndarray) -> np.ndarray:
    """Per-column (cell type) min-max normalization."""
    mins = X.min(axis=0)
    maxs = X.max(axis=0)
    rng = maxs - mins
    rng[rng == 0] = 1.0
    return (X - mins) / rng


def _batch_correct_quantile(comp: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Per-batch quantile normalization to match global distribution."""
    n, p = comp.shape
    corrected = comp.copy()
    for j in range(p):
        global_sorted = np.sort(comp[:, j])
        global_quantiles = np.linspace(0, 1, n)
        for b in np.unique(labels):
            mask = labels == b
            if mask.sum() < 3:
                continue
            batch_vals = comp[mask, j]
            ranks = np.argsort(np.argsort(batch_vals)).astype(float)
            ranks = ranks / max(len(batch_vals) - 1, 1)
            corrected[mask, j] = np.interp(ranks, global_quantiles, global_sorted)
    corrected = np.maximum(corrected, 0.0)
    row_sums = corrected.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    return corrected / row_sums


def _batch_correct_rank(comp: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Rank-based inverse normal mapping: rank within batch, map to global."""
    n, p = comp.shape
    global_sorted = {j: np.sort(comp[:, j]) for j in range(p)}
    corrected = comp.copy()
    for j in range(p):
        for b in np.unique(labels):
            mask = labels == b
            if mask.sum() < 3:
                continue
            batch_vals = comp[mask, j]
            n_b = len(batch_vals)
            ranks = np.argsort(np.argsort(batch_vals)).astype(float)
            quantiles = (ranks + 0.5) / n_b
            global_indices = np.clip((quantiles * n).astype(int), 0, n - 1)
            corrected[mask, j] = global_sorted[j][global_indices]
    corrected = np.maximum(corrected, 0.0)
    row_sums = corrected.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    return corrected / row_sums


def _apply_batch_correction(freq_dict: dict, batch_labels: dict,
                            method: str = "quantile") -> dict:
    """Apply batch correction to composition data.

    Parameters
    ----------
    freq_dict : dict of {group_name: patients x cell_types array}
    batch_labels : dict of {group_name: int array of batch labels}
    method : "quantile", "rank", or "none"

    Returns
    -------
    Corrected freq_dict (same structure).
    """
    if method == "none":
        return freq_dict

    correct_fn = _batch_correct_quantile if method == "quantile" else _batch_correct_rank
    corrected = {}
    for g, comp in freq_dict.items():
        if g in batch_labels and len(np.unique(batch_labels[g])) > 1:
            corrected[g] = correct_fn(comp, batch_labels[g])
        else:
            corrected[g] = comp
    return corrected
