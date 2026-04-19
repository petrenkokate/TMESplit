"""Numerical constants and combiner registry."""

EPSILON = 1e-10

VALID_COMBINERS = ("cauchy", "bonferroni", "fisher_permutation", "fisher_asymptotic")

_COMBINER_DISPLAY = {
    "cauchy": "Cauchy (ACAT)",
    "bonferroni": "Bonferroni",
    "fisher_permutation": "Fisher (permutation)",
    "fisher_asymptotic": "Fisher (asymptotic)",
}
