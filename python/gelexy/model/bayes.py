import gc

import numpy as np
from formulae import design_matrices

from gelexy import BayesAlphabet
from gelexy._core import _BayesModel
from gelexy.data import read_fam
from gelexy.utils import intersection

from ..formula import format_formula
from .base import ModelMaker, get_fixed_levels


class BayesModel(_BayesModel):
    pass


class make_bayes(ModelMaker):
    def make(self, formula: str, bfile: str, bayes_type: BayesAlphabet):
        """
        Create a Linear Mixed Model from the specified formula and genetic relationship matrices.

        :param formula: A formula string specifying the model. The left-hand side specifies the response variable,
                       and the right-hand side specifies the fixed effects. Example: "y ~ x1 + x2"
        :param bfile: the prefix of plink bed files
        :returns: A GBLUP object containing the response vector, design matrix, genetic relationship
                  matrices, and random effect names.
        :rtype: GBLUP

        :raises ValueError: If the fixed effects contain missing values.
        """

        formula = format_formula(formula)
        response = formula.split("~")[0].strip()

        data = self.data.copy()
        data = self._dropna(data, response)

        # intersect with genotype data
        gid = read_fam(bfile, self._iid_only)
        data, pid, gid = intersection(data, gid, self._iid_only)
        self.logger.info("")

        design_matrix = design_matrices(formula, data, na_action="error")

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

        model._add_genetic_effect(
            self._iid_only,
            False,
            bfile,
            pid,
            gid,
            bayes_type,
        )

        gc.collect()

        return model
