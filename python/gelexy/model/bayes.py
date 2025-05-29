import gc

import numpy as np
import pandas as pd
from formulae import design_matrices

from gelexy import BayesAlphabet
from gelexy._core import BayesModel, _BedReader, _sp_dense_dot

from .base import ModelMakerBase, create_design_matrix_genetic
from .formula_parser import FormulaParser


class make_bayes(ModelMakerBase):
    def make(
        self, formula: str, genotypes: dict[str], bayes_type: BayesAlphabet
    ):
        """
        Create a Linear Mixed Model from the specified formula and genetic relationship matrices.

        Parameters
        ----------
        formula : str
            A formula string specifying the model. The left-hand side specifies the response variable,
            and the right-hand side specifies the fixed effects. Example: "y ~ x1 + x2"
        grm : dict[str, pd.DataFrame | str | Path]
            A dictionary mapping random effect names to their corresponding genetic relationship matrices.
            The matrices can be provided as DataFrames or paths to files containing the matrices.

        Returns
        -------
        GBLUP
            A GBLUP object containing the response vector, design matrix, genetic relationship
            matrices, and random effect names.

        Raises
        ------
        ValueError
            If the fixed effects contain missing values.
        """
        fparser = FormulaParser(formula, list(genotypes.keys()))
        data = self.data
        data = self._clear_data(data, fparser.response)

        for name, genotype in genotypes.items():
            genotypes[name] = load_genotype(genotype)

        data, design_matrix_genetic, dropped_ids, genotypes = (
            create_design_matrix_genetic(data, genotypes, symmetric=False)
        )

        design_mat = design_matrices(fparser.common, data, na_action="error")

        model = BayesModel(
            fparser.format_common,
            np.asfortranarray(design_mat.response),
        )

        model.add_fixed_effect(
            list(design_mat.common.terms),
            self._get_fixed_levels(design_mat.common.terms),
            np.asfortranarray(design_mat.common.design_matrix),
        )

        if design_mat.group is not None:
            for name, matrix in design_mat.group.terms.items():
                model.add_random_effect(
                    "(" + name + ")", np.asfortranarray(matrix.data)
                )

        for term in fparser.genetic_terms:
            genotype = genotypes[term.genetic].to_numpy()
            genotype = _sp_dense_dot(design_matrix_genetic, genotype)
            assert genotype.flags["F_CONTIGUOUS"], (
                "Genotype matrix must be Fortran contiguous."
            )
            model.add_genetic_effect(
                term.name,
                genotype,
                bayes_type,
            )

        # model._dropped_ids = dropped_ids
        # model._formula = formula
        gc.collect()
        return model


def load_genotype(bed_file: str):
    reader = _BedReader(bed_file, int(1e9))
    return pd.DataFrame(
        reader.read_chunk(True),
        index=reader.individuals,
        columns=reader.snps,
    )
