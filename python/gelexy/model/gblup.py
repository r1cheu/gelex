import gc
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from formulae import design_matrices
from scipy.sparse import csc_matrix

from gelexy._core import _GBLUP
from gelexy.data.grm import load_grm

from .base import ModelMakerBase, create_design_matrix_genetic
from .formula_parser import FormulaParser


class GBLUP(_GBLUP):
    @property
    def U(self):
        return pd.DataFrame(
            self._U, index=self._individuals, columns=self.random_effect_names
        )

    def save(self, path: str | Path | None = None):
        """
        Save the model to a file.

        Parameters
        ----------
        path : str | Path | None
            Path to the file where the model will be saved.
        """
        if path is None:
            path = (
                f"{self._lhs}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.model"
            )

        with h5py.File(path, "w") as f:
            f.create_dataset("beta", data=self.beta)
            f.create_dataset("sigma", data=self.sigma)
            f.create_dataset("proj_y", data=self._proj_y)
            f.create_dataset("dropped_ids", data=self._dropped_ids)
            f.attrs["formula"] = self.formula()


class make_gblup(ModelMakerBase):
    def make(self, formula: str, grms: dict) -> GBLUP:
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
        fparser = FormulaParser(formula, list(grms.keys()))
        data = self.data
        data = self._clear_data(data, fparser.response)
        grms = self._load_grms(grms)

        data, design_matrix_genetic, dropped_ids, grms = (
            create_design_matrix_genetic(data, grms)
        )

        design_mat = design_matrices(fparser.common, data, na_action="error")

        model = GBLUP(
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
                    "(" + name + ")", csc_matrix(matrix.data)
                )

        for term in fparser.genetic_terms:
            model.add_genetic_effect(
                term.name,
                design_matrix_genetic,
                np.asfortranarray(grms[term.genetic]),
            )

        for term in fparser.gxe_terms:
            dm = design_matrices("0+" + term.env, data, na_action="error")
            model.add_gxe_effect(
                term.name,
                design_matrix_genetic,
                np.asfortranarray(grms[term.genetic]),
                np.asfortranarray(dm.common),
            )
        model.add_residual()

        model._dropped_ids = dropped_ids
        model._formula = formula
        gc.collect()
        return model

    @staticmethod
    def _load_grms(grms):
        for name, grm in grms.items():
            if isinstance(grm, str | Path):
                grms[name] = load_grm(grm)
            elif isinstance(grm, pd.DataFrame):
                grms[name] = grm
            else:
                msg = f"GRM for {name} must be a path to a file or a DataFrame."
                raise ValueError(msg)
        return grms
