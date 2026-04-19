"""scvi-tools-style ``TMESplit`` class.

Wraps the functional ``ecosplit`` + ``test_significance`` API around an
AnnData that holds patient x cell-type compositions.
"""

from __future__ import annotations

from collections.abc import Sequence

import anndata as ad
import numpy as np
import pandas as pd

from tmesplit._core import ecosplit
from tmesplit._testing import test_significance

_UNS_KEY = "_tmesplit_config"


class TMESplit:
    """Estimator for group-specific multicellular coordination programs.

    Follows the scvi-tools pattern: register fields on an AnnData via
    :meth:`setup_anndata`, construct with the AnnData, then call
    :meth:`train` to fit NMF and run the hierarchical permutation test.

    Expected AnnData shape
    ----------------------
    ``adata`` is a patient-level composition object:
      * ``adata.n_obs`` = number of patients
      * ``adata.n_vars`` = number of cell types
      * ``adata.X`` = proportions (rows should sum to 1, but the fit
        applies per-group min-max normalization so exact unit sums are
        not required)
      * ``adata.obs`` has one categorical column identifying the group
        (registered via :meth:`setup_anndata`)

    For cell-level data, aggregate first via
    :func:`tmesplit.pp.aggregate_to_composition`.

    Parameters
    ----------
    adata
        Patient-level composition AnnData. Must have gone through
        :meth:`setup_anndata`.
    **ecosplit_kwargs
        Forwarded to :func:`tmesplit.ecosplit` at ``train`` time
        (e.g. ``k_hint``, ``n_runs``, ``beta_loss``, ``match_threshold``).

    Examples
    --------
    >>> TMESplit.setup_anndata(adata, group_key="condition")
    >>> model = TMESplit(adata, beta_loss=0.5)
    >>> model.train(n_perm=2000, random_state=42)
    >>> model.p_value
    0.003
    """

    def __init__(self, adata: ad.AnnData, **ecosplit_kwargs):
        cfg = adata.uns.get(_UNS_KEY)
        if cfg is None or "group_key" not in cfg:
            raise RuntimeError(
                "AnnData has not been registered. Call "
                "TMESplit.setup_anndata(adata, group_key=...) first."
            )
        self.adata = adata
        self._group_key: str = cfg["group_key"]
        self._ecosplit_kwargs = dict(ecosplit_kwargs)
        self._fitted = False
        self._result: dict | None = None
        self._perm: dict | None = None

    @classmethod
    def setup_anndata(cls, adata: ad.AnnData, group_key: str) -> None:
        """Register the obs column that identifies biological groups.

        Mirrors ``scvi.model.SCVI.setup_anndata``: stores the registration
        in ``adata.uns[_tmesplit_config]`` so subsequent constructors can
        find it without requiring the user to pass ``group_key`` again.
        """
        if group_key not in adata.obs.columns:
            raise KeyError(f"{group_key!r} not found in adata.obs")
        n_groups = adata.obs[group_key].astype(str).nunique()
        if n_groups < 2:
            raise ValueError(
                f"Need at least 2 distinct groups in adata.obs[{group_key!r}], "
                f"got {n_groups}"
            )
        adata.uns[_UNS_KEY] = {"group_key": group_key, "n_groups": int(n_groups)}

    @classmethod
    def from_matrix(
        cls,
        X: np.ndarray,
        groups: Sequence,
        cell_type_names: Sequence[str] | None = None,
        sample_names: Sequence[str] | None = None,
        group_name: str = "group",
        **ecosplit_kwargs,
    ) -> TMESplit:
        """Build a ``TMESplit`` from a plain matrix + group labels.

        Used by the R basilisk bridge and by Python users without AnnData.
        Internally constructs an AnnData wrapper and calls
        :meth:`setup_anndata`.
        """
        X_arr = np.asarray(X, dtype=np.float64)
        if X_arr.ndim != 2:
            raise ValueError("X must be 2-dimensional (patients x cell types)")
        groups_arr = np.asarray(groups)
        if groups_arr.shape[0] != X_arr.shape[0]:
            raise ValueError(
                f"groups has {groups_arr.shape[0]} entries but X has "
                f"{X_arr.shape[0]} rows"
            )

        var_index = (
            list(cell_type_names)
            if cell_type_names is not None
            else [f"ct_{i}" for i in range(X_arr.shape[1])]
        )
        obs_index = (
            list(sample_names)
            if sample_names is not None
            else [f"sample_{i}" for i in range(X_arr.shape[0])]
        )

        obs = pd.DataFrame({group_name: groups_arr.astype(str)}, index=obs_index)
        var = pd.DataFrame(index=var_index)

        adata = ad.AnnData(X=X_arr, obs=obs, var=var)
        cls.setup_anndata(adata, group_key=group_name)
        return cls(adata, **ecosplit_kwargs)

    # ----- fit & predict -----

    def _freq_dict(self) -> dict:
        labels = self.adata.obs[self._group_key].astype(str).to_numpy()
        X = self.adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float64)
        unique = sorted(np.unique(labels).tolist())
        return {g: X[labels == g] for g in unique}

    def train(
        self,
        n_perm: int = 2000,
        n_jobs: int = 1,
        random_state: int = 42,
        verbose: bool = False,
        fit_kwargs: dict | None = None,
        test_kwargs: dict | None = None,
    ) -> TMESplit:
        """Fit NMF and run the 3-level hierarchical permutation test.

        Parameters
        ----------
        n_perm
            Number of permutations for the hierarchical test. Default 2000.
        n_jobs
            Parallel jobs for the permutation loop (joblib). ``-1`` = all cores,
            ``1`` = serial (default, byte-reproducible).
        """
        freq_dict = self._freq_dict()

        fk = {"seed": random_state, "verbose": verbose, **self._ecosplit_kwargs,
              **(fit_kwargs or {})}
        self._result = ecosplit(freq_dict, **fk)

        tk = {"seed": random_state, "verbose": verbose, "n_perm": n_perm,
              "n_jobs": n_jobs, **(test_kwargs or {})}
        self._perm = test_significance(freq_dict, result=self._result, **tk)

        self._fitted = True
        return self

    def _require_fitted(self):
        if not self._fitted or self._result is None or self._perm is None:
            raise RuntimeError("Call .train(...) before accessing results.")

    # ----- accessors (attribute style, scvi/sklearn-ish) -----

    @property
    def groups(self) -> list[str]:
        self._require_fitted()
        return self._result["groups"]

    @property
    def k_shared(self) -> int:
        self._require_fitted()
        return self._result["k_shared"]

    @property
    def k_specific(self) -> dict:
        self._require_fitted()
        return self._result["k_specific"]

    @property
    def W_shared(self) -> np.ndarray:
        self._require_fitted()
        return self._result["W_shared"]

    @property
    def W_specific(self) -> dict:
        self._require_fitted()
        return self._result["W_specific"]

    @property
    def H_fractions(self) -> dict:
        self._require_fitted()
        return self._result["H_fractions"]

    @property
    def p_value(self) -> float:
        self._require_fitted()
        return self._perm["p_value"]

    @property
    def p_value_per_group(self) -> dict:
        self._require_fitted()
        return self._perm["p_value_per_group"]

    @property
    def programs(self) -> pd.DataFrame:
        """Program-level table with per-program p-values and FDR."""
        self._require_fitted()
        table = self._perm["programs_table"]
        if not table:
            return pd.DataFrame()
        return pd.DataFrame(table)

    @property
    def result_(self) -> dict:
        """Raw ``ecosplit()`` result dict."""
        self._require_fitted()
        return self._result

    @property
    def permutation_(self) -> dict:
        """Raw ``test_significance()`` result dict."""
        self._require_fitted()
        return self._perm

    def __repr__(self) -> str:
        try:
            state = "fitted" if self._fitted else "unfitted"
            n_obs, n_vars = self.adata.shape
            extras = f", k_shared={self.k_shared}, p={self.p_value:.4g}" if self._fitted else ""
            return (f"TMESplit({state}, n_obs={n_obs}, n_vars={n_vars}, "
                    f"group_key={self._group_key!r}{extras})")
        except Exception:
            return object.__repr__(self)
