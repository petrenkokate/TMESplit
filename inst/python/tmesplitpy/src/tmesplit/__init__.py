"""TMESplit -- formal testing of group-specific multicellular coordination programs.

Two layers of API:

Functional (stable, drives the regression fixtures):
  ecosplit()                   Fit per-group NMF + match + classify programs.
  test_significance()          3-level hierarchical permutation test.
  minmax_normalize()           Per-column min-max normalization helper.
  evaluate_on_simulation()     Benchmarking helper.

Class (scvi-tools-style, AnnData-first):
  TMESplit                     Estimator class with setup_anndata / train / accessors.
  TMESplit.from_matrix(...)    Build from a plain patients x cell-types matrix.

Submodules:
  tmesplit.pp                  Preprocessing helpers (e.g. aggregate_to_composition).
  tmesplit.pl                  Plotting helpers (matplotlib; requires [plot] extra).
"""

from tmesplit import pl, pp
from tmesplit._batch import minmax_normalize
from tmesplit._core import ecosplit
from tmesplit._evaluation import evaluate_on_simulation
from tmesplit._testing import test_significance
from tmesplit._version import __version__
from tmesplit.model import TMESplit

__all__ = [
    "TMESplit",
    "ecosplit",
    "test_significance",
    "minmax_normalize",
    "evaluate_on_simulation",
    "pp",
    "pl",
    "__version__",
]
