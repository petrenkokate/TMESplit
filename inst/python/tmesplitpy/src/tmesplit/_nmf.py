"""NMF engine, parallel-analysis rank selection, NNLS projection, and VE helpers."""

from __future__ import annotations

import warnings

import numpy as np
from scipy.linalg import svd
from scipy.optimize import nnls
from sklearn.decomposition import NMF

from tmesplit._constants import EPSILON


def _run_nmf(V, k, n_runs=10, max_iter=500, seed=42, beta_loss=0.5):
    """V: features x samples. Returns W (p, k), H (k, n), error.

    beta_loss: 2 = Frobenius, 1 = KL divergence, 0.5 = beta divergence.
    KL (beta=1) and beta=0.5 give better program recovery on compositional data.
    """
    V_pos = np.maximum(V, EPSILON)
    best_err = np.inf
    best_W = best_H = None
    rng = np.random.default_rng(seed)
    solver = 'mu' if beta_loss != 2 else 'cd'
    init = 'nndsvda' if beta_loss == 2 else 'nndsvda'
    for _run in range(n_runs):
        rs = int(rng.integers(0, 2**31))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = NMF(n_components=k, init=init, max_iter=max_iter,
                    solver=solver, beta_loss=beta_loss, random_state=rs)
            Wsk = m.fit_transform(V_pos.T)
            Hsk = m.components_
            err = m.reconstruction_err_
        if err < best_err:
            best_err = err
            best_W = Hsk.T   # (p, k)
            best_H = Wsk.T   # (k, n)
    return best_W, best_H, best_err


def _parallel_analysis(X: np.ndarray, n_perm: int = 100, percentile: float = 90,
                       seed: int = 42) -> int:
    """Select number of components via parallel analysis on eigenvalues.

    Compares eigenvalues of the data to those of column-permuted data.
    Returns the number of eigenvalues exceeding the permutation threshold.

    X: patients x features (already normalized).
    """
    rng = np.random.default_rng(seed)
    n, p = X.shape
    X_c = X - X.mean(axis=0)
    s_real = svd(X_c, compute_uv=False)
    eig_real = s_real ** 2 / (n - 1)

    eig_perm = np.zeros((n_perm, min(n, p)))
    for i in range(n_perm):
        X_perm = X.copy()
        for col in range(p):
            rng.shuffle(X_perm[:, col])
        s_p = svd(X_perm - X_perm.mean(axis=0), compute_uv=False)
        eig_perm[i, :len(s_p)] = s_p ** 2 / (n - 1)

    eig_thresh = np.percentile(eig_perm, percentile, axis=0)
    k = int(np.sum(eig_real[:len(eig_thresh)] > eig_thresh[:len(eig_real)]))
    return max(k, 2)


def _select_rank(X: np.ndarray, n_runs: int = 5, max_iter: int = 500,
                 recon_threshold: float = 0.04, pa_percentile: float = 90,
                 seed: int = 42) -> int:
    """Two-step rank selection: parallel analysis + reconstruction improvement.

    1. PA gives base k (conservative).
    2. If NMF reconstruction improves by >threshold when going from k to k+1,
       bump to k+1.

    Validated: 82.8% exact, 97.8% within +/-1 on 9 scenarios x 10 seeds.
    Atlas: 100% exact (Female=6, Male=5).
    """
    k_base = _parallel_analysis(X, percentile=pa_percentile, seed=seed)

    _, _, err_k = _run_nmf(X.T, k_base, n_runs=n_runs, max_iter=max_iter, seed=seed)
    _, _, err_k1 = _run_nmf(X.T, k_base + 1, n_runs=n_runs, max_iter=max_iter, seed=seed)

    improvement = (err_k - err_k1) / (err_k + EPSILON)
    if improvement > recon_threshold:
        return k_base + 1
    return k_base


def _project_onto_basis(V, H_basis):
    """Project V (n x p) onto H_basis (k x p) via NNLS. Returns W (n x k)."""
    n = V.shape[0]
    k = H_basis.shape[0]
    W = np.zeros((n, k))
    for i in range(n):
        W[i], _ = nnls(H_basis.T, V[i])
    return W


def _best_cross_corr(w, W_other):
    best = 0.0
    for j in range(W_other.shape[1]):
        c = np.corrcoef(w, W_other[:, j])[0, 1]
        if not np.isnan(c):
            best = max(best, abs(c))
    return best


def _program_ve(W, H, V):
    k = W.shape[1]
    total = np.sum(V ** 2)
    ve = np.zeros(k)
    for j in range(k):
        ve[j] = np.sum(np.outer(W[:, j], H[j]) ** 2) / (total + EPSILON)
    return ve
