"""Three-level hierarchical permutation test (Global / Per-group / Per-program)."""

from __future__ import annotations

import numpy as np
from joblib import Parallel, delayed
from scipy.stats import chi2

from tmesplit._batch import _apply_batch_correction, minmax_normalize
from tmesplit._constants import _COMBINER_DISPLAY, EPSILON, VALID_COMBINERS
from tmesplit._core import ecosplit
from tmesplit._matching import _match_programs_multigroup
from tmesplit._nmf import (
    _best_cross_corr,
    _program_ve,
    _project_onto_basis,
    _run_nmf,
    _select_rank,
)

try:
    from numba import njit
    _HAS_NUMBA = True
except ImportError:  # pragma: no cover
    _HAS_NUMBA = False

    def njit(*args, **kwargs):
        def _wrap(fn):
            return fn
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return _wrap


@njit(cache=True, fastmath=False)
def _rank_columns(T_all: np.ndarray) -> np.ndarray:
    """For each column j and row ii, compute (1 + #{T_all[:,j] >= T_all[ii,j]}) / (1+n).

    Hot kernel for Fisher-under-permutation combination: O(n^2 * ncols) but
    numba-friendly (tight loops, numeric arrays). Called once per test_significance
    invocation with n = n_perm + 1.
    """
    n, k = T_all.shape
    out = np.empty((n, k), dtype=np.float64)
    denom = 1.0 + n
    for j in range(k):
        for ii in range(n):
            thr = T_all[ii, j]
            c = 0
            for r in range(n):
                if T_all[r, j] >= thr:
                    c += 1
            out[ii, j] = (1.0 + c) / denom
    return out


def _compute_test_statistics(V_dict, groups, k_per_group, n_runs=3,
                              max_iter=200, beta_loss=0.5, seed=42):
    """Run per-group NMF and compute test statistics for specificity.

    V_dict: {group: patients x features (already per-group normalized)}
    Supports G >= 2 groups.

    Test statistics:
      recon_asymmetry: within-group recon improvement minus cross-group
      ve_ratio: variance explained by specific factor in own vs other group
      projection_gap: Cohen's d of specific factor scores between groups

    Returns dict of {stat_name: value}.
    """
    p = V_dict[groups[0]].shape[1]

    W, H = {}, {}
    for g in groups:
        Vt = V_dict[g].T
        W[g], H[g], _ = _run_nmf(Vt, k_per_group[g], n_runs=n_runs,
                                   max_iter=max_iter, seed=seed, beta_loss=beta_loss)

    shared_comps, corrs, um = _match_programs_multigroup(W, groups, threshold=0.5)
    k_shared = len(shared_comps)

    specific = {}
    for g in groups:
        spec = []
        ve = _program_ve(W[g], H[g], V_dict[g].T)
        for i in um.get(g, []):
            max_br = 0.0
            for g2 in groups:
                if g2 != g:
                    max_br = max(max_br, _best_cross_corr(W[g][:, i], W[g2]))
            if max_br < 0.4 and ve[i] >= 0.02:
                spec.append(i)
        specific[g] = spec

    n_specific = {g: len(specific[g]) for g in groups}
    has_specific = sum(n_specific.values()) > 0

    if not has_specific:
        return {
            'has_specific': False, 'n_specific': n_specific, 'k_shared': k_shared,
            'recon_asymmetry': 0.0, 've_ratio': 1.0, 'projection_gap': 0.0,
            'programs': [],
            'recon_asymmetry_per': np.array([], dtype=float),
            've_ratio_per': np.array([], dtype=float),
            'projection_gap_per': np.array([], dtype=float),
        }

    recon_asymmetries = []
    for g in groups:
        for idx in specific[g]:
            w_spec = W[g][:, idx]
            h_spec_basis = w_spec.reshape(1, -1)

            shared_cols = [comp[g] for comp in shared_comps if g in comp]
            if shared_cols:
                H_shared = W[g][:, shared_cols].T
                H_aug = np.vstack([H_shared, h_spec_basis])
            else:
                H_shared = np.zeros((0, p))
                H_aug = h_spec_basis

            if H_shared.shape[0] > 0:
                W_sh = _project_onto_basis(V_dict[g], H_shared)
                mse_shared_own = np.mean((V_dict[g] - W_sh @ H_shared) ** 2)
            else:
                mse_shared_own = np.mean(V_dict[g] ** 2)
            W_aug = _project_onto_basis(V_dict[g], H_aug)
            mse_aug_own = np.mean((V_dict[g] - W_aug @ H_aug) ** 2)
            delta_own = mse_shared_own - mse_aug_own

            delta_cross_list = []
            for g2 in groups:
                if g2 == g:
                    continue
                shared_cols_g2 = [comp[g2] for comp in shared_comps if g2 in comp]
                if shared_cols_g2:
                    H_sh_g2 = W[g2][:, shared_cols_g2].T
                    H_aug_g2 = np.vstack([H_sh_g2, h_spec_basis])
                else:
                    H_sh_g2 = np.zeros((0, p))
                    H_aug_g2 = h_spec_basis
                if H_sh_g2.shape[0] > 0:
                    W_sh2 = _project_onto_basis(V_dict[g2], H_sh_g2)
                    mse_sh_cross = np.mean((V_dict[g2] - W_sh2 @ H_sh_g2) ** 2)
                else:
                    mse_sh_cross = np.mean(V_dict[g2] ** 2)
                W_aug2 = _project_onto_basis(V_dict[g2], H_aug_g2)
                mse_aug_cross = np.mean((V_dict[g2] - W_aug2 @ H_aug_g2) ** 2)
                delta_cross_list.append(mse_sh_cross - mse_aug_cross)
            recon_asymmetries.append(delta_own - np.mean(delta_cross_list))

    recon_asymmetry = np.max(recon_asymmetries) if recon_asymmetries else 0.0

    ve_ratios = []
    for g in groups:
        for idx in specific[g]:
            w_spec = W[g][:, idx]
            var_own = np.var(V_dict[g] @ w_spec)
            vars_cross = [np.var(V_dict[g2] @ w_spec) for g2 in groups if g2 != g]
            ve_ratios.append(var_own / (np.mean(vars_cross) + EPSILON))
    ve_ratio = np.max(ve_ratios) if ve_ratios else 1.0

    projection_gaps = []
    for g in groups:
        for idx in specific[g]:
            w_spec = W[g][:, idx]
            w_norm = w_spec / (np.linalg.norm(w_spec) + EPSILON)
            scores_own = V_dict[g] @ w_norm
            scores_cross = np.concatenate([V_dict[g2] @ w_norm
                                            for g2 in groups if g2 != g])
            pooled_std = np.sqrt((np.var(scores_own) + np.var(scores_cross)) / 2 + EPSILON)
            projection_gaps.append(abs(np.mean(scores_own) - np.mean(scores_cross)) / pooled_std)
    projection_gap = np.max(projection_gaps) if projection_gaps else 0.0

    program_ids = [(g, int(idx)) for g in groups for idx in specific[g]]

    return {
        'has_specific': True, 'n_specific': n_specific, 'k_shared': k_shared,
        'recon_asymmetry': float(recon_asymmetry),
        've_ratio': float(ve_ratio),
        'projection_gap': float(projection_gap),
        'programs': program_ids,
        'recon_asymmetry_per': np.asarray(recon_asymmetries, dtype=float),
        've_ratio_per': np.asarray(ve_ratios, dtype=float),
        'projection_gap_per': np.asarray(projection_gaps, dtype=float),
    }


def _one_perm(V_perm, groups, k_per_group_obs, rerank_in_perm, ps,
              n_runs, max_iter, recon_threshold, beta_loss, stats_names, _floor):
    """Single-permutation worker (module-level for joblib picklability)."""
    try:
        if rerank_in_perm:
            k_per_group_perm = {}
            for g in groups:
                k_per_group_perm[g] = _select_rank(
                    V_perm[g], n_runs=n_runs, max_iter=max_iter,
                    recon_threshold=recon_threshold, seed=ps)
        else:
            k_per_group_perm = k_per_group_obs
        stats = _compute_test_statistics(V_perm, groups, k_per_group_perm,
                                          n_runs=n_runs, max_iter=max_iter,
                                          beta_loss=beta_loss, seed=ps)
    except Exception:
        return None

    out = {s: stats[s] for s in stats_names}
    out_per_group = {g: {s: _floor[s] for s in stats_names} for g in groups}
    pool = {s: [] for s in stats_names}

    if stats['has_specific']:
        perm_prog_groups = [pg for (pg, _) in stats['programs']]
        for g in groups:
            mask = np.array([pg == g for pg in perm_prog_groups], dtype=bool)
            if mask.any():
                for s in stats_names:
                    out_per_group[g][s] = float(np.max(stats[f'{s}_per'][mask]))
        for s in stats_names:
            pool[s] = stats[f'{s}_per'].tolist()

    return {'stats': out, 'per_group': out_per_group, 'pool': pool,
            'has_specific': stats['has_specific']}


def test_significance(freq_dict, result=None, n_perm=500, n_runs=3,
                      max_iter=200, beta_loss=0.5,
                      batch_labels=None, batch_correction="none",
                      combiners="cauchy",
                      recon_threshold=0.04,
                      rerank_in_perm=True,
                      n_jobs=1,
                      seed=42, verbose=True):
    """Permutation test for group-specific programs (hierarchical, 3-level).

    Shuffles group labels, re-runs per-group normalization + NMF + matching,
    and computes test statistics. Per-group normalization is re-applied inside
    each permutation to preserve label exchangeability.

    Three test statistics (one per geometric aspect of group specificity):
      recon_asymmetry : within-group reconstruction improvement minus cross-group
        (shape -- does the factor improve fit for its own group more than others?)
      ve_ratio        : variance explained by the specific factor in its own group
        relative to other groups (dispersion)
      projection_gap  : Cohen's d of specific factor scores between groups (location)

    Combiners (choose via `combiners` argument):
      cauchy             -- Cauchy / ACAT (Liu 2019 AJHG). DEFAULT.
      bonferroni         -- min(1, 3 * min(marginal p-values)). Conservative.
      fisher_permutation -- Fisher's statistic ranked in permutation null.
      fisher_asymptotic  -- Asymptotic chi2_6 Fisher. ANTI-CONSERVATIVE.

    When batch_labels are provided, uses Freedman-Lane residual permutation.

    Hierarchical output (closed-testing discipline):
      Level 1 -- Global: `p_value` = primary combiner applied to max-stats
      Level 2 -- Per-group: `p_value_per_group[g]` = primary combiner on group-
                stratified max-stats. Only interpret if Level 1 is significant.
      Level 3 -- Per-program: `programs_table` gives DESeq-like rows with
                primary combiner applied per detected specific program.
    """
    combiner_list = [combiners] if isinstance(combiners, str) else list(combiners)
    if not combiner_list:
        raise ValueError("combiners must be non-empty")
    for c in combiner_list:
        if c not in VALID_COMBINERS:
            raise ValueError(
                f"Unknown combiner: {c!r}. Valid options: {VALID_COMBINERS}"
            )
    primary_combiner = combiner_list[0]

    groups = sorted(freq_dict.keys())
    G = len(groups)
    assert G >= 2, f"Need at least 2 groups, got {G}"

    freq_corrected = freq_dict
    if batch_labels is not None and batch_correction != "none":
        freq_corrected = _apply_batch_correction(freq_dict, batch_labels, batch_correction)

    if result is None:
        result = ecosplit(freq_corrected, n_runs=max(n_runs * 2, 5), seed=seed,
                          verbose=verbose, beta_loss=beta_loss)

    k_per_group = result['k_per_group']

    if verbose:
        print(f"\nPermutation test (k_per_group={k_per_group})")

    V_obs = {g: minmax_normalize(freq_corrected[g]) for g in groups}

    if verbose:
        print("Computing observed test statistics...")
    T_obs = _compute_test_statistics(V_obs, groups, k_per_group,
                                      n_runs=max(n_runs * 2, 5),
                                      max_iter=500, beta_loss=beta_loss, seed=seed)

    if verbose:
        if T_obs['has_specific']:
            print(f"  recon_asymmetry = {T_obs['recon_asymmetry']:.6f}")
            print(f"  ve_ratio        = {T_obs['ve_ratio']:.3f}")
            print(f"  projection_gap  = {T_obs['projection_gap']:.3f}")
        else:
            print("  No specific programs found")

    if not T_obs['has_specific']:
        if verbose:
            print("  Skipping permutation -- no specific programs to test")
        p_values_empty = {'recon_asymmetry': 1.0, 've_ratio': 1.0, 'projection_gap': 1.0}
        for c in combiner_list:
            p_values_empty[c] = 1.0
        return {
            'p_value': 1.0,
            'p_value_per_group': {g: 1.0 for g in groups},
            'p_values_per_group_details': {g: {} for g in groups},
            'p_values_all': p_values_empty,
            'programs_table': [],
            'T_obs': T_obs, 'n_perm': 0,
            'combiners': combiner_list,
            'primary_combiner': primary_combiner,
        }

    all_patients = np.vstack([freq_dict[g] for g in groups])
    n_per_group = {g: freq_dict[g].shape[0] for g in groups}
    n_total = all_patients.shape[0]
    split_indices = []
    cumsum = 0
    for g in groups:
        cumsum += n_per_group[g]
        split_indices.append(cumsum)

    all_batch_labels = None
    if batch_labels is not None:
        all_batch_labels = np.concatenate([batch_labels[g] for g in groups])

    rng = np.random.default_rng(seed + 7777)
    stats_names = ['recon_asymmetry', 've_ratio', 'projection_gap']
    _floor = {'recon_asymmetry': 0.0, 've_ratio': 1.0, 'projection_gap': 0.0}

    perm_orders = []
    perm_seeds = []
    for _ in range(n_perm):
        perm_orders.append(rng.permutation(n_total))
        perm_seeds.append(int(rng.integers(0, 2**31)))

    def _prepare_V_perm(perm):
        patients_perm = all_patients[perm]
        slices = np.split(patients_perm, split_indices[:-1])
        if all_batch_labels is not None and batch_correction != "none":
            batch_perm = all_batch_labels[perm]
            batch_slices = np.split(batch_perm, split_indices[:-1])
            perm_dict = {g: slices[j] for j, g in enumerate(groups)}
            batch_perm_dict = {g: batch_slices[j] for j, g in enumerate(groups)}
            perm_corrected = _apply_batch_correction(perm_dict, batch_perm_dict, batch_correction)
            return {g: minmax_normalize(perm_corrected[g]) for g in groups}
        return {g: minmax_normalize(slices[j]) for j, g in enumerate(groups)}

    V_perms = [_prepare_V_perm(p) for p in perm_orders]

    if verbose:
        print(f"Running {n_perm} permutations (n_jobs={n_jobs})...")

    results = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(_one_perm)(
            V_perms[i], groups, k_per_group, rerank_in_perm, perm_seeds[i],
            n_runs, max_iter, recon_threshold, beta_loss, stats_names, _floor,
        )
        for i in range(n_perm)
    )

    T_perm = {s: [] for s in stats_names}
    T_perm_per_group = {g: {s: [] for s in stats_names} for g in groups}
    pooled_null = {s: [] for s in stats_names}

    for r in results:
        if r is None:
            continue
        for s in stats_names:
            T_perm[s].append(r['stats'][s])
        for g in groups:
            for s in stats_names:
                T_perm_per_group[g][s].append(r['per_group'][g][s])
        if r['has_specific']:
            for s in stats_names:
                pooled_null[s].extend(r['pool'][s])

    n_done = len(T_perm['recon_asymmetry'])
    if verbose:
        p_ra = (1 + sum(t >= T_obs['recon_asymmetry'] for t in T_perm['recon_asymmetry'])) / (n_done + 1)
        p_ve = (1 + sum(t >= T_obs['ve_ratio'] for t in T_perm['ve_ratio'])) / (n_done + 1)
        print(f"  {n_done}/{n_perm} succeeded: p_recon={p_ra:.4f}, p_ve={p_ve:.4f}")

    def _combine_pvalues(T_obs_vec, T_perm_arrays):
        """Compute marginal p-values + all 4 combiners for a single test."""
        marg = {}
        for s in stats_names:
            arr = np.asarray(T_perm_arrays[s])
            count = int(np.sum(arr >= T_obs_vec[s]))
            marg[s] = (1 + count) / (1 + len(arr))

        T_stack = np.column_stack([np.asarray(T_perm_arrays[s], dtype=np.float64)
                                    for s in stats_names])
        T_all = np.vstack([T_stack, np.array([[T_obs_vec[s] for s in stats_names]],
                                              dtype=np.float64)])
        n_all = T_all.shape[0]
        p_rank = _rank_columns(T_all)
        fisher_all = -2.0 * np.sum(np.log(np.maximum(p_rank, 1e-10)), axis=1)
        p_fperm = float(
            (1 + int(np.sum(fisher_all[:-1] >= fisher_all[-1]))) / (1 + (n_all - 1))
        )

        fstat_asym = float(np.sum(
            -2 * np.log(np.maximum([marg[s] for s in stats_names], 1e-10))
        ))
        p_asym = float(1 - chi2.cdf(fstat_asym, df=6))

        p_clip = np.clip([marg[s] for s in stats_names], 1e-10, 1 - 1e-10)
        cstat = float(np.mean(np.tan((0.5 - p_clip) * np.pi)))
        p_cauchy = float(0.5 - np.arctan(cstat) / np.pi)

        p_bonf = float(min(1.0, 3.0 * min(marg.values())))

        combiners_out = {
            "cauchy": p_cauchy,
            "bonferroni": p_bonf,
            "fisher_permutation": p_fperm,
            "fisher_asymptotic": p_asym,
        }
        return marg, combiners_out

    T_obs_L1 = {s: T_obs[s] for s in stats_names}
    p_marg, combiners_L1 = _combine_pvalues(T_obs_L1, T_perm)

    p_values = dict(p_marg)
    for c in combiner_list:
        p_values[c] = combiners_L1[c]

    p_value_per_group = {}
    p_values_per_group_details = {}
    obs_programs = T_obs['programs']
    obs_program_groups = [pg for (pg, _) in obs_programs]
    for g in groups:
        mask_obs = np.array([pg == g for pg in obs_program_groups], dtype=bool)
        if not mask_obs.any():
            p_value_per_group[g] = 1.0
            details = {s: 1.0 for s in stats_names}
            for c in combiner_list:
                details[c] = 1.0
            p_values_per_group_details[g] = details
            continue
        T_obs_g = {
            s: float(np.max(T_obs[f'{s}_per'][mask_obs])) for s in stats_names
        }
        marg_g, combiners_g = _combine_pvalues(T_obs_g, T_perm_per_group[g])
        p_value_per_group[g] = combiners_g[primary_combiner]
        details = dict(marg_g)
        for c in combiner_list:
            details[c] = combiners_g[c]
        p_values_per_group_details[g] = details

    def _bh_fdr(pvals):
        pvals = np.asarray(pvals, dtype=float)
        n = len(pvals)
        if n == 0:
            return pvals
        order = np.argsort(pvals)
        ranked = pvals[order]
        q = ranked * n / (np.arange(n) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        q = np.clip(q, 0, 1)
        out = np.empty(n)
        out[order] = q
        return out

    level3_combiner = primary_combiner
    level3_fallback_note = None
    if level3_combiner == "fisher_permutation":
        level3_combiner = "cauchy"
        level3_fallback_note = (
            "Fisher-permutation is not available at Level 3 (no per-program "
            "permutation null); falling back to Cauchy (ACAT)."
        )
    if level3_combiner == "fisher_asymptotic":
        pass

    def _level3_combine(marginal_p_list, combiner_name):
        p_arr = np.clip(list(marginal_p_list), 1e-10, 1 - 1e-10)
        if combiner_name == "cauchy":
            cstat = float(np.mean(np.tan((0.5 - p_arr) * np.pi)))
            return float(0.5 - np.arctan(cstat) / np.pi)
        elif combiner_name == "bonferroni":
            return float(min(1.0, 3.0 * min(p_arr)))
        elif combiner_name == "fisher_asymptotic":
            fstat = float(np.sum(-2 * np.log(p_arr)))
            return float(1 - chi2.cdf(fstat, df=6))
        else:
            raise ValueError(f"Unsupported Level 3 combiner: {combiner_name}")

    programs_table = []
    for i, (g, idx) in enumerate(obs_programs):
        row = {'group': g, 'program_idx': int(idx)}
        marginal_p = {}
        for s in stats_names:
            obs_val = float(T_obs[f'{s}_per'][i])
            pool = np.asarray(pooled_null[s])
            n_null = len(pool)
            if n_null == 0:
                p = 1.0
            else:
                count = int(np.sum(pool >= obs_val))
                p = (1 + count) / (1 + n_null)
            row[s] = obs_val
            row[f'{s}_p'] = float(p)
            marginal_p[s] = p
        row['combined_p'] = _level3_combine(
            [marginal_p[s] for s in stats_names], level3_combiner
        )
        programs_table.append(row)

    for g in groups:
        g_indices = [i for i, r in enumerate(programs_table) if r['group'] == g]
        if not g_indices:
            continue
        ps_g = np.array([programs_table[i]['combined_p'] for i in g_indices])
        qs_g = _bh_fdr(ps_g)
        for i, q in zip(g_indices, qs_g, strict=False):
            programs_table[i]['fdr'] = float(q)

    primary_p = p_values[primary_combiner]

    def _sig_marker(p):
        return ("***" if p < 0.001 else
                "**"  if p < 0.01  else
                "*"   if p < 0.05  else "")

    if verbose:
        print(f"\n{'='*60}")
        print(f"HIERARCHICAL P-VALUES ({n_done} permutations)")
        print(f"Primary combiner: {_COMBINER_DISPLAY[primary_combiner]}")
        print(f"{'='*60}")
        print("LEVEL 1 -- Global")
        for s in stats_names:
            p = p_values[s]
            print(f"  {s:<28s} T={T_obs[s]:.4f}  p={p:.4f} {_sig_marker(p)}")
        for c in combiner_list:
            p = p_values[c]
            marker = "  [primary]" if c == primary_combiner else ""
            label = _COMBINER_DISPLAY[c]
            print(f"  {label:<28s}          p={p:.4f} {_sig_marker(p)}{marker}")

        print(f"\nLEVEL 2 -- Per-group ({_COMBINER_DISPLAY[primary_combiner]})")
        if primary_p >= 0.05:
            print("  (Level 1 not significant -- closed-testing: do not interpret)")
        for g in groups:
            p = p_value_per_group[g]
            print(f"  {g:<28s}          p={p:.4f} {_sig_marker(p)}")

        l3_label = _COMBINER_DISPLAY[level3_combiner]
        print(f"\nLEVEL 3 -- Per-program ({l3_label} + BH-FDR within group)")
        if level3_fallback_note:
            print(f"  note: {level3_fallback_note}")
        if primary_p >= 0.05:
            print("  (Level 1 not significant -- closed-testing: do not interpret)")
        for row in programs_table:
            gp = row['group']
            drill_ok = p_value_per_group[gp] < 0.05
            mark = "" if drill_ok else "  (group n.s.)"
            cp, fdr = row['combined_p'], row.get('fdr', 1.0)
            print(f"  {gp}/prog{row['program_idx']:<4d}  "
                  f"ve={row['ve_ratio']:.2f}  ve_p={row['ve_ratio_p']:.4f}  "
                  f"combined_p={cp:.4f}  FDR={fdr:.4f} {_sig_marker(fdr)}{mark}")

    return {
        'p_value': primary_p,
        'p_value_per_group': p_value_per_group,
        'p_values_per_group_details': p_values_per_group_details,
        'p_values_all': p_values,
        'programs_table': programs_table,
        'T_obs': T_obs,
        'n_perm': n_done,
        'combiners': combiner_list,
        'primary_combiner': primary_combiner,
        'level3_combiner': level3_combiner,
    }
