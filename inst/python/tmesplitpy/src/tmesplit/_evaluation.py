"""Simulation evaluation helpers (program recovery, assignment accuracy)."""

from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment

from tmesplit._constants import EPSILON


def _soft_assignment_accuracy(H_true, H_est):
    """Correlation between true and estimated H fractions, Hungarian-matched."""
    H_t = H_true / (np.sum(H_true, axis=1, keepdims=True) + EPSILON)
    H_e = H_est / (np.sum(H_est, axis=1, keepdims=True) + EPSILON)

    kt, ke = H_t.shape[1], H_e.shape[1]
    corr = np.zeros((kt, ke))
    for i in range(kt):
        for j in range(ke):
            c = np.corrcoef(H_t[:, i], H_e[:, j])[0, 1]
            corr[i, j] = c if not np.isnan(c) else 0
    nm = max(kt, ke)
    cost = np.zeros((nm, nm))
    cost[:kt, :ke] = -corr
    ri, ci = linear_sum_assignment(cost)
    matched = [corr[r, c] for r, c in zip(ri, ci, strict=False) if r < kt and c < ke]
    return float(np.mean(matched)) if matched else 0.0


def _program_recovery(Wt, We):
    if Wt.shape[1] == 0 or We.shape[1] == 0:
        return 0.0
    kt, ke = Wt.shape[1], We.shape[1]
    corr = np.zeros((kt, ke))
    for i in range(kt):
        for j in range(ke):
            c = np.corrcoef(Wt[:, i], We[:, j])[0, 1]
            corr[i, j] = abs(c) if not np.isnan(c) else 0.0
    n = max(kt, ke)
    cost = np.zeros((n, n))
    cost[:kt, :ke] = -corr
    ri, ci = linear_sum_assignment(cost)
    m = [corr[r, c] for r, c in zip(ri, ci, strict=False) if r < kt and c < ke]
    return float(np.mean(m)) if m else 0.0


def _assignment_accuracy(ta, ea):
    if len(ta) != len(ea):
        return 0.0
    n = len(ta)
    tl, el = np.unique(ta), np.unique(ea)
    nm = max(len(tl), len(el))
    cost = np.zeros((nm, nm))
    tm = {l: i for i, l in enumerate(tl)}
    em = {l: i for i, l in enumerate(el)}
    for i in range(n):
        ti, ei = tm.get(ta[i]), em.get(ea[i])
        if ti is not None and ei is not None:
            cost[ti, ei] -= 1
    ri, ci = linear_sum_assignment(cost)
    return float(-cost[ri, ci].sum() / n)


def evaluate_on_simulation(sim_data, result):
    kst = sim_data['W_shared'].shape[1]
    groups = result['groups']
    kse = result['k_shared']

    sr = _program_recovery(sim_data['W_shared'], result['W_shared'])

    spec_r, kspt = {}, {}
    for gi, g in enumerate(groups):
        kspt[g] = sim_data['W_specific'][gi].shape[1]
        wt, we = sim_data['W_specific'][gi], result['W_specific'][g]
        if wt.shape[1] > 0 and we.shape[1] > 0:
            spec_r[g] = _program_recovery(wt, we)
        elif wt.shape[1] == 0 and we.shape[1] == 0:
            spec_r[g] = 1.0
        else:
            spec_r[g] = 0.0
    ms = np.mean(list(spec_r.values())) if spec_r else 1.0

    aa = {}
    for gi, g in enumerate(groups):
        aa[g] = _assignment_accuracy(
            sim_data['assignments'][gi], result['patient_assignments'][g])
    ma = np.mean(list(aa.values()))

    sa = {}
    for gi, g in enumerate(groups):
        H_true = sim_data['H'][gi]
        H_est = result['H_fractions'][g]
        sa[g] = _soft_assignment_accuracy(H_true, H_est)
    msa = np.mean(list(sa.values()))

    return {
        'shared_recovery': sr, 'specific_recovery': spec_r,
        'mean_specific_recovery': ms,
        'assignment_accuracy': aa, 'mean_assignment_accuracy': ma,
        'soft_assignment_accuracy': sa, 'mean_soft_assignment': msa,
        'k_shared_true': kst, 'k_shared_est': kse,
        'k_shared_correct': int(kse == kst),
        'k_specific_true': kspt,
        'k_specific_est': result['k_specific'],
        'k_specific_correct': {g: int(result['k_specific'][g] == kspt[g]) for g in groups},
        'composite': (sr + ms + msa) / 3.0,
        'composite_soft': (sr + ms + msa) / 3.0,
        'composite_hard': (sr + ms + ma) / 3.0,
    }
