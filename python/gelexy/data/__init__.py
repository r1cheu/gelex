"""This module initializes the dataset module, provide tools for handle Genotype, Phenotype and Calculation of Genetic Relation Matrix."""

from .grm import load_grm, make_grm
from .reader import read_fam, read_phenotype

__all__ = [
    "load_grm",
    "make_grm",
    "read_fam",
    "read_phenotype",
]
