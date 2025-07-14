import gc
from pathlib import Path

import numpy as np
import pandas as pd
from formulae import design_matrices

from gelexy import BayesAlphabet
from gelexy._core import _BayesModel, _BedReader, _sp_dense_dot
from gelexy.utils.aligner import Matcher

from ..formula import Formula
from .base import ModelMakerBase, get_fixed_levels, make_design_matrix


class BayesModel(_BayesModel):
    pass


class make_bayes(ModelMakerBase):
    def make(
        self, formula: str, genotypes: dict[str], bayes_type: BayesAlphabet
    ):
        """
        Create a Linear Mixed Model from the specified formula and genetic relationship matrices.

        :param formula: A formula string specifying the model. The left-hand side specifies the response variable,
                       and the right-hand side specifies the fixed effects. Example: "y ~ x1 + x2"
        :type formula: str
        :param genotypes: A dictionary mapping random effect names to their corresponding genetic relationship matrices.
                          The matrices can be provided as DataFrames or paths to files containing the matrices.
        :type genotypes: dict[str, pd.DataFrame | str | Path]

        :returns: A GBLUP object containing the response vector, design matrix, genetic relationship
                  matrices, and random effect names.
        :rtype: GBLUP

        :raises ValueError: If the fixed effects contain missing values.
        """
        fparser = Formula(formula, list(genotypes.keys()))

        data = self.data.copy()
        data = self._dropna(data, fparser.response)
        genotypes = self._load_genotype(genotypes)

        self.matcher = Matcher(data, genotypes)
        data = self.matcher.phenotype
        genotypes = self.matcher.genotypes

        design_matrix = design_matrices(fparser.common, data, na_action="error")
        design_mat_genotype = make_design_matrix(
            data["id"],
            self.matcher.common_order,
        )

        model = BayesModel(
            fparser.formula,
            np.asfortranarray(design_matrix.response),
        )

        model._add_fixed_effect(
            list(design_matrix.common.terms),
            get_fixed_levels(design_matrix.common.terms),
            np.asfortranarray(design_matrix.common.design_matrix),
        )

        if design_matrix.group is not None:
            for name, matrix in design_matrix.group.terms.items():
                model._add_random_effect(
                    "(" + name + ")", np.asfortranarray(matrix)
                )

        for term in fparser.genetic_terms:
            genotype = np.asfortranarray(genotypes[term.genetic])
            genotype = _sp_dense_dot(design_mat_genotype, genotype)
            model._add_genetic_effect(
                term.name,
                genotype,
                bayes_type,
            )
        gc.collect()
        return model

    def _load_genotype(self, genotypes: dict[str, pd.DataFrame | str | Path]):
        """
        Load genotype data from the provided paths or DataFrames.

        :param genotypes: A dictionary mapping random effect names to their corresponding genetic relationship matrices.
                          The matrices can be provided as DataFrames or paths to files containing the matrices.
        :type genotypes: dict[str, pd.DataFrame | str | Path]

        :returns: A dictionary of loaded genotype DataFrames.
        :rtype: dict[str, pd.DataFrame]
        """
        for name, genotype in genotypes.items():
            if isinstance(genotype, str | Path):
                genotypes[name] = load_genotype(genotype)
            elif isinstance(genotype, pd.DataFrame):
                genotypes[name] = genotype
            else:
                msg = f"genotype for {name} must be a path to a plink bed or a DataFrame."
                raise ValueError(msg)
        return genotypes


def load_genotype(bed_file: str):
    reader = _BedReader(bed_file, int(1e9))
    genotype = pd.DataFrame(
        reader.read_chunk(),
        index=reader.individuals,
        columns=reader.snps,
    )
    genotype.index.name = "id"
    return genotype
