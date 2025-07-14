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

from ..formula import Formula
from .base import ModelMakerBase, get_fixed_levels, make_design_matrix


class GBLUP(_GBLUP):
    @property
    def U(self):
        return pd.DataFrame(
            self._U, index=self._individuals, columns=self.random_effect_names
        )

    @property
    def train_sample(self):
        if hasattr(self, "_train_sample"):
            return self._train_sample
        msg = "Train sample is not set. Are you sure this is a trained model?"
        raise ValueError(msg)

    @train_sample.setter
    def train_sample(self, value):
        self._train_sample = value

    def save(self, path: str | Path | None = None):
        """
        Save the model parameters to an HDF5 file.
        :param path: The file path to save the model. If None, a default name is generated.
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
        fparser = Formula(formula, list(grms.keys()))

        data = self.data.copy()
        data = self._dropna(data, fparser.response)
        grms = self._load_grms(grms)

        self.matcher = Matcher(data, grms, axis=[0, 1])
        data = self.matcher.phenotype
        grms = self.matcher.genotypes

        design_matrix = design_matrices(fparser.common, data, na_action="error")
        design_mat_grm = make_design_matrix(
            data["id"],
            self.matcher.common_order,
        )

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

        model.train_sample = deepcopy(self.matcher.common_order)
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
