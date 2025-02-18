"""This module initializes the phenx package, which is designed to provide tools and utilities for genomics selection."""

# import pandas as pd

# from .dataset import Genotypes, Phenotypes, grm, intersect
# from .model import LinearMixedModel

# pd.options.mode.copy_on_write = True

# __all__ = ["Genotypes", "Phenotypes", "grm", "LinearMixedModel", "intersect"]

from ._chenx import Estimator
from .data import load_grm, make_grm, read_table
from .model import make_model

__all__ = ["Estimator", "load_grm", "make_grm", "make_model", "read_table"]
