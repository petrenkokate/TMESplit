"""Main ecosplit() function: per-group NMF + cross-group matching + program selection."""

from __future__ import annotations

from itertools import product as iterproduct

import numpy as np

from tmesplit._batch import _apply_batch_correction, minmax_normalize
from tmesplit._constants import EPSILON
from tmesplit._matching import _match_programs_multigroup
from tmesplit._nmf import _best_cross_corr, _program_ve, _run_nmf, _select_rank


def ecosplit(freq_dict: dict[str, np.ndarray],
             cell_type_names: list[str] | None = None,
             k_hint: int | None = None,
             n_runs: int = 10,
             match_threshold: float = 0.5,
             min_ve_specific: float = 0.02,
             recon_threshold: float = 0.04,
             beta_loss: float = 0.5,
             batch_labels: dict[str, np.ndarray] | None = None,
             batch_correction: str = "none",
             seed: int = 42,
             verbose: bool = True) -> dict:
    """Run EcoSplit v2 with support for G >= 2 groups.

    Parameters
    ----------
    freq_dict : dict of {group_name: np.ndarray (patients x cell_types)}
    k_hint : int, optional
        If provided, use as per-group NMF rank (try k_hint and k_hint+1).
        If None, automatic rank selection via parallel analysis + reconstruction.
    n_runs : int
        NMF restarts. Default 10.
    match_threshold : float
        Min |r| for shared programs. Default 0.5.
    min_ve_specific : float
        Min variance explained for a specific program. Default 0.02.
    recon_threshold : float
        Min reconstruction improvement (fraction) to justify k+1 over k.
        Only used in automatic mode. Default 0.04.
    batch_labels : dict of {group_name: int array}, optional
        Per-group batch/dataset labels. If provided, batch correction is
        applied before NMF.
    batch_correction : str
        Batch correction method: "quantile", "rank", or "none" (default).
    """
    groups = sorted(freq_dict.keys())
    G = len(groups)
    assert G >= 2, f"Need at least 2 groups, got {G}"

    p = freq_dict[groups[0]].shape[1]
    if cell_type_names is None:
        cell_type_names = [f'ct_{i}' for i in range(p)]

    if batch_labels is not None and batch_correction != "none":
        if verbose:
            print(f"Applying batch correction: {batch_correction}")
        freq_dict = _apply_batch_correction(freq_dict, batch_labels, batch_correction)

    V = {}
    V_norm = {}
    for g in groups:
        V_norm[g] = minmax_normalize(freq_dict[g])
        V[g] = V_norm[g].T

    if k_hint is not None:
        k_per_group = {g: k_hint for g in groups}
        k_try = {g: [k_hint, k_hint + 1] for g in groups}
        if verbose:
            print(f"Using k_hint={k_hint}, testing k and k+1")
    else:
        k_per_group = {}
        k_try = {}
        for g in groups:
            k_auto = _select_rank(V_norm[g], n_runs=n_runs,
                                  recon_threshold=recon_threshold, seed=seed)
            k_per_group[g] = k_auto
            k_try[g] = [max(k_auto - 1, 2), k_auto]
            if verbose:
                print(f"  Auto rank {g}: k={k_auto}")

    W_cache = {g: {} for g in groups}
    H_cache = {g: {} for g in groups}

    for g in groups:
        for k in k_try[g]:
            W_g, H_g, err = _run_nmf(V[g], k, n_runs=n_runs, seed=seed,
                                      beta_loss=beta_loss)
            W_cache[g][k] = W_g
            H_cache[g][k] = H_g

    best_combo = None
    best_score = -1.0

    k_lists = [k_try[g] for g in groups]
    for combo in iterproduct(*k_lists):
        k_combo = {g: k for g, k in zip(groups, combo, strict=False)}

        W_combo = {g: W_cache[g][k_combo[g]] for g in groups}
        shared_comps, corrs, um = _match_programs_multigroup(
            W_combo, groups, threshold=match_threshold)

        n_shared = len(shared_comps)
        mean_r = float(np.mean(corrs)) if corrs else 0.0

        n_spec = {}
        for g in groups:
            n_spec_g = 0
            ve = _program_ve(W_combo[g], H_cache[g][k_combo[g]], V[g])
            other_groups = [og for og in groups if og != g]
            for i in um[g]:
                max_br = 0.0
                for og in other_groups:
                    br = _best_cross_corr(W_combo[g][:, i], W_combo[og])
                    max_br = max(max_br, br)
                if max_br < match_threshold - 0.1 and ve[i] >= min_ve_specific:
                    n_spec_g += 1
            n_spec[g] = n_spec_g

        n_noise = sum(k_combo[g] - n_shared - n_spec[g] for g in groups)
        score = n_shared * mean_r - 0.5 * n_noise

        if verbose:
            k_str = ','.join(str(k_combo[g]) for g in groups)
            sp_str = '+'.join(str(n_spec[g]) for g in groups)
            print(f"  k=({k_str}): {n_shared} shared (r={mean_r:.3f}), "
                  f"spec={sp_str}, noise={n_noise}, score={score:.3f}")

        if score > best_score:
            best_score = score
            best_combo = k_combo

    if verbose:
        k_str = ','.join(str(best_combo[g]) for g in groups)
        print(f"\n  Best: k=({k_str})")

    W_final = {g: W_cache[g][best_combo[g]] for g in groups}
    H_final = {g: H_cache[g][best_combo[g]] for g in groups}

    shared_comps, corrs, unmatched = _match_programs_multigroup(
        W_final, groups, threshold=match_threshold)

    specific_final = {}
    for g in groups:
        sig = []
        ve = _program_ve(W_final[g], H_final[g], V[g])
        other_groups = [og for og in groups if og != g]
        for i in unmatched[g]:
            max_br = 0.0
            for og in other_groups:
                br = _best_cross_corr(W_final[g][:, i], W_final[og])
                max_br = max(max_br, br)
            if max_br < match_threshold - 0.1 and ve[i] >= min_ve_specific:
                sig.append(i)
                if verbose:
                    print(f"  specific[{g}]: prog[{i}] VE={ve[i]:.4f}, best_r={max_br:.3f}")
        specific_final[g] = sig

    k_shared = len(shared_comps)
    k_specific = {g: len(specific_final[g]) for g in groups}

    if verbose:
        print(f"\n  k_shared={k_shared}, k_specific={k_specific}")
        for comp, r in zip(shared_comps, corrs, strict=False):
            parts = ' <-> '.join(f"{g}[{comp[g]}]" for g in groups)
            print(f"    shared: {parts} r={r:.3f}")

    W_shared = np.zeros((p, k_shared))
    sidx = {}
    for s, comp in enumerate(shared_comps):
        ws = []
        w_ref = W_final[groups[0]][:, comp[groups[0]]]
        for g in groups:
            w_g = W_final[g][:, comp[g]]
            if np.corrcoef(w_ref, w_g)[0, 1] < 0:
                w_g = -w_g
            ws.append(w_g)
        W_shared[:, s] = np.mean(ws, axis=0)
        for g in groups:
            sidx[(g, comp[g])] = s

    W_spec_out, spidx = {}, {}
    for g in groups:
        ks = k_specific[g]
        if ks > 0:
            W_spec_out[g] = np.column_stack([W_final[g][:, i] for i in specific_final[g]])
            for s, i in enumerate(specific_final[g]):
                spidx[(g, i)] = s
        else:
            W_spec_out[g] = np.zeros((p, 0))

    H_sh, H_sp = {}, {}
    for g in groups:
        ng = freq_dict[g].shape[0]
        kg = W_final[g].shape[1]
        hs = np.zeros((k_shared, ng))
        for i in range(kg):
            if (g, i) in sidx:
                hs[sidx[(g, i)]] = H_final[g][i]
        H_sh[g] = hs
        ksg = k_specific[g]
        hp = np.zeros((ksg, ng))
        for i in range(kg):
            if (g, i) in spidx:
                hp[spidx[(g, i)]] = H_final[g][i]
        H_sp[g] = hp

    pa, H_frac = {}, {}
    for g in groups:
        hf = np.vstack([H_sh[g], H_sp[g]])
        pa[g] = np.argmax(hf, axis=0)
        row_sums = hf.sum(axis=0, keepdims=True) + EPSILON
        H_frac[g] = (hf / row_sums).T
        if verbose:
            ng = freq_dict[g].shape[0]
            u, c = np.unique(pa[g], return_counts=True)
            print(f"\n  '{g}': {len(u)} ecotypes")
            for ui, ci in zip(u, c, strict=False):
                l = f"shared_{ui}" if ui < k_shared else f"specific_{g}_{ui-k_shared}"
                print(f"    {l}: {ci} ({100*ci/ng:.1f}%)")

    return {
        'k_shared': k_shared, 'k_specific': k_specific,
        'k_per_group': best_combo,
        'W_shared': W_shared, 'W_specific': W_spec_out,
        'W_per_group': W_final,
        'H_shared': H_sh, 'H_specific': H_sp,
        'H_fractions': H_frac,
        'patient_assignments': pa,
        'cell_type_names': cell_type_names, 'groups': groups,
        'shared_matches': shared_comps,
    }
