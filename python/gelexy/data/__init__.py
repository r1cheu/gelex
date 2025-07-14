"""This module initializes the dataset module, provide tools for handle Genotype, Phenotype and Calculation of Genetic Relation Matrix."""

from .grm import load_genotype, load_grm, make_grm
from .reader import read_pheno

__all__ = ["load_genotype", "load_grm", "make_grm", "read_pheno"]
