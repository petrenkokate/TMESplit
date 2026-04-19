"""Cross-correlation matching of NMF programs across groups."""

from __future__ import annotations

from collections import defaultdict

import numpy as np
from scipy.optimize import linear_sum_assignment


def _match_programs(W_a, W_b, threshold=0.5):
    """Match programs between two groups using |correlation|."""
    ka, kb = W_a.shape[1], W_b.shape[1]
    if ka == 0 or kb == 0:
        return [], list(range(ka)), list(range(kb))
    corr = np.zeros((ka, kb))
    for i in range(ka):
        for j in range(kb):
            c = np.corrcoef(W_a[:, i], W_b[:, j])[0, 1]
            corr[i, j] = c if not np.isnan(c) else 0.0
    n = max(ka, kb)
    cost = np.zeros((n, n))
    cost[:ka, :kb] = -np.abs(corr)
    ri, ci = linear_sum_assignment(cost)
    shared, ma, mb = [], set(), set()
    for r, c in zip(ri, ci, strict=False):
        if r < ka and c < kb and abs(corr[r, c]) >= threshold:
            shared.append((r, c, float(corr[r, c])))
            ma.add(r); mb.add(c)
    return shared, [i for i in range(ka) if i not in ma], [j for j in range(kb) if j not in mb]


def _match_programs_multigroup(W_dict, groups, threshold=0.5):
    """Match programs across G>=2 groups via all-pairs pairwise matching.

    A program is "shared" if it participates in a connected component that
    spans ALL groups in the pairwise match graph. Otherwise it's unmatched.

    Returns
    -------
    shared_components : list of dict {group: program_index}
        Each entry maps every group to its matched program index.
    mean_corrs : list of float
        Mean |correlation| for each shared component.
    unmatched : dict of {group: list of unmatched program indices}
    """
    G = len(groups)
    if G == 2:
        g0, g1 = groups
        shared_pairs, um0, um1 = _match_programs(W_dict[g0], W_dict[g1], threshold)
        components = [{g0: ia, g1: ib} for ia, ib, _ in shared_pairs]
        corrs = [abs(r) for _, _, r in shared_pairs]
        return components, corrs, {g0: um0, g1: um1}

    pair_matches = {}
    for i in range(G):
        for j in range(i + 1, G):
            gi, gj = groups[i], groups[j]
            shared, _, _ = _match_programs(W_dict[gi], W_dict[gj], threshold)
            pair_matches[(gi, gj)] = shared

    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for g in groups:
        for k_idx in range(W_dict[g].shape[1]):
            parent[(g, k_idx)] = (g, k_idx)

    for (gi, gj), matches in pair_matches.items():
        for ia, ib, _r in matches:
            union((gi, ia), (gj, ib))

    components_raw = defaultdict(set)
    for g in groups:
        for k_idx in range(W_dict[g].shape[1]):
            root = find((g, k_idx))
            components_raw[root].add((g, k_idx))

    shared_components = []
    mean_corrs = []
    matched_nodes = set()

    for root, nodes in components_raw.items():
        groups_in_comp = set(g for g, _ in nodes)
        if len(groups_in_comp) == G and len(nodes) == G:
            comp = {g: idx for g, idx in nodes}
            corr_vals = []
            for i in range(G):
                for j in range(i + 1, G):
                    gi, gj = groups[i], groups[j]
                    c = np.corrcoef(W_dict[gi][:, comp[gi]],
                                    W_dict[gj][:, comp[gj]])[0, 1]
                    if not np.isnan(c):
                        corr_vals.append(abs(c))
                    else:
                        corr_vals.append(0.0)
            shared_components.append(comp)
            mean_corrs.append(float(np.mean(corr_vals)))
            matched_nodes.update(nodes)

    unmatched = {}
    for g in groups:
        kg = W_dict[g].shape[1]
        unmatched[g] = [i for i in range(kg) if (g, i) not in matched_nodes]

    return shared_components, mean_corrs, unmatched
