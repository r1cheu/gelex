import gc
from pathlib import Path

import numpy as np
import pandas as pd
from formulae import design_matrices
from scipy.sparse import csc_matrix

from gelexy._core import _GBLUP
from gelexy.data.grm import load_grm
from gelexy.utils import align_gblup

from ..formula import Formula
from .base import ModelMakerBase, get_fixed_levels, make_design_matrix


class GBLUP(_GBLUP):
    @property
    def U(
        self,
    ):
        if not hasattr(self, "_genotype_id"):
            msg = ""
            raise AttributeError(msg)

        return pd.DataFrame(self.u, index=self._genotype_id)


class make_gblup(ModelMakerBase):
    def make(self, formula: str, grms: dict) -> GBLUP:
        fparser = Formula(formula, list(grms.keys()))

        data = self.data.copy()
        grms = self._load_grms(grms)
        data, genotype_id = align_gblup(data, grms)

        data = self._dropna(data, fparser.response)

        design_matrix = design_matrices(fparser.common, data, na_action="error")
        design_mat_grm = make_design_matrix(data["id"], genotype_id)

        model = GBLUP(
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

        model._formula = formula
        model._genotype_id = genotype_id
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
