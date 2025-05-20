"""This module initializes the gelexy package, which is designed to provide tools and utilities for genomics selection."""

# import pandas as pd

# from .dataset import Genotypes, Phenotypes, grm, intersect
# from .model import GBLUP

# pd.options.mode.copy_on_write = True

# __all__ = ["Genotypes", "Phenotypes", "grm", "GBLUP", "intersect"]

from ._core import MCMC, BayesAlphabet, Estimator
from .data import load_grm, make_grm, read_table
from .model import make_bayes, make_gblup

__all__ = [
    "MCMC",
    "BayesAlphabet",
    "Estimator",
    "load_grm",
    "make_bayes",
    "make_gblup",
    "make_grm",
    "read_table",
]
