"""Shared helpers for the ``pl`` submodule."""

from __future__ import annotations


def _require_matplotlib():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
        import seaborn as sns  # noqa: F401
    except ImportError as e:
        raise ImportError(
            "Plotting requires matplotlib + seaborn. "
            "Install with `pip install tmesplitpy[plot]`."
        ) from e


def _require_fitted(model):
    if not getattr(model, "_fitted", False):
        raise RuntimeError("Call `model.train(...)` before plotting.")
