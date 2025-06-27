import gc
from copy import deepcopy
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from formulae import design_matrices
from scipy.sparse import csc_matrix

from gelexy._core import _GBLUP
from gelexy.data.grm import load_grm
from gelexy.utils.aligner import Matcher

from .base import (
    ModelMakerBase,
    make_design_matrix,
)
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

        data = self.data.copy()
        data = self._dropna(data, fparser.response)
        grms = self._load_grms(grms)

        self.matcher = Matcher(data, grms, axis=[0, 1])
        data = self.matcher.phenotype
        grms = self.matcher.genotypes

        design_mat = design_matrices(fparser.common, data, na_action="error")
        design_mat_grm = make_design_matrix(
            data["id"],
            self.matcher.common_order,
        )

        model = GBLUP(
            fparser.formula,
            np.asfortranarray(design_mat.response),
        )

        model._add_fixed_effect(
            list(design_mat.common.terms),
            self._get_fixed_levels(design_mat.common.terms),
            np.asfortranarray(design_mat.common.design_matrix),
        )

        if design_mat.group is not None:
            for name, matrix in design_mat.group.terms.items():
                model._add_random_effect(
                    "(" + name + ")", csc_matrix(matrix.data)
                )

        for term in fparser.genetic_terms:
            model._add_genetic_effect(
                term.name,
                design_mat_grm,
                np.asfortranarray(grms[term.genetic]),
            )

        for term in fparser.gxe_terms:
            dm = design_matrices("0+" + term.env, data, na_action="error")
            model._add_gxe_effect(
                term.name,
                design_mat_grm,
                np.asfortranarray(grms[term.genetic]),
                np.asfortranarray(dm.common),
            )

        model._add_residual()

        model._genotype_order = deepcopy(self.matcher.common_order)
        model._formula = formula
        gc.collect()
        return model

    def _load_grms(self, grms):
        for name, grm in grms.items():
            if isinstance(grm, str | Path):
                grms[name] = load_grm(grm)
            elif isinstance(grm, pd.DataFrame):
                grms[name] = grm
            else:
                msg = f"GRM for {name} must be a path to a file or a DataFrame."
                raise ValueError(msg)
        return grms
