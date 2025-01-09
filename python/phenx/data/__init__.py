"""This module initializes the dataset module, provide tools for handle Genotype, Phenotype and Calculation of Genetic Relation Matrix."""

from .grm import _load_grm, make_grm

__all__ = ["_load_grm", "make_grm"]
