"""This module initializes the gelexy package, which is designed to provide tools and utilities for genomics selection."""

# import pandas as pd

# from .dataset import Genotypes, Phenotypes, grm, intersect
# from .model import GBLUP

# pd.options.mode.copy_on_write = True

# __all__ = ["Genotypes", "Phenotypes", "grm", "GBLUP", "intersect"]

from ._core import MCMC, BayesAlphabet, Estimator, MCMCParams
from .data import load_genotype, load_grm, make_grm, read_pheno
from .model import GBLUP, BayesModel, make_bayes, make_gblup
from .predictor import BayesPredictor

__all__ = [
    "GBLUP",
    "MCMC",
    "BayesAlphabet",
    "BayesModel",
    "BayesPredictor",
    "Estimator",
    "MCMCParams",
    "load_genotype",
    "load_grm",
    "make_bayes",
    "make_gblup",
    "make_grm",
    "read_pheno",
]
