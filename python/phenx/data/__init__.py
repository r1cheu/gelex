"""This module initializes the dataset module, provide tools for handle Genotype, Phenotype and Calculation of Genetic Relation Matrix."""

from .grm import load_grm, make_grm

__all__ = ["load_grm", "make_grm"]
