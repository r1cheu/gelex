"""This module initializes the phenx package, which is designed to provide tools and utilities for genomics selection."""

from .dataset import Genotypes, Phenotypes, grm

# from .optim import REML

__all__ = ["Genotypes", "Phenotypes", "grm"]
