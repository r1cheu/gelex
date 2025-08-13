import gc

import numpy as np
from formulae import design_matrices

from gelexy import BayesAlphabet, load_genotype
from gelexy._core import _BayesModel, _sp_dense_dot
from gelexy.utils import align_bayes

from ..formula import format_formula
from .base import ModelMakerBase, get_fixed_levels, make_design_matrix


class BayesModel(_BayesModel):
    pass


class make_bayes(ModelMakerBase):
    def make(self, formula: str, bfile: str, bayes_type: BayesAlphabet):
        """
        Create a Linear Mixed Model from the specified formula and genetic relationship matrices.

        :param formula: A formula string specifying the model. The left-hand side specifies the response variable,
                       and the right-hand side specifies the fixed effects. Example: "y ~ x1 + x2"
        :type formula: str
        :param genotypes: A dictionary mapping random effect names to their corresponding genetic relationship matrices.
                          The matrices can be provided as DataFrames or paths to files containing the matrices.

        :returns: A GBLUP object containing the response vector, design matrix, genetic relationship
                  matrices, and random effect names.
        :rtype: GBLUP

        :raises ValueError: If the fixed effects contain missing values.
        """

        formula = format_formula(formula)
        response = formula.split("~")[0].strip()

        data = self.data.copy()
        data = self._dropna(data, response)
        genotype = load_genotype(bfile)

        data, genotype = align_bayes(data, genotype)

        design_matrix = design_matrices(formula, data, na_action="error")
        design_mat_genotype = make_design_matrix(
            data["id"],
            genotype.index.to_list(),
        )

        model = BayesModel(
            formula,
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

        genotype = np.asfortranarray(genotype.to_numpy())
        genotype = _sp_dense_dot(design_mat_genotype, genotype)
        model._add_genetic_effect(
            "a",
            genotype,
            bayes_type,
        )
        gc.collect()
        return model
