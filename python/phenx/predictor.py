from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from formulaic import Formula

from phenx.utils.path import valid_path

from ._chenx import _Predictor
from .data import read_table
from .model import LinearMixedModel, LinearMixedModelParams
from .model.io import load_params


class Predictor(_Predictor):
    def predict(
        self, test_bed: str | Path, data: str | Path | pd.DataFrame | None = None
    ):
        test_bed = valid_path(test_bed, (".bed",))
        u = self._compute_u(test_bed)

        covariates = self._set_covariates(self.test_individuals, data)
        pred = u.sum(axis=1, keepdims=True) + self._compute_covariates(covariates)
        return pd.DataFrame(
            np.hstack([pred, u]), columns=[self._params.lhs, *self._random_effect_names]
        )

    def _set_covariates(
        self, test_individuals: list[str], data: str | Path | pd.DataFrame | None = None
    ):
        rhs = self._params.rhs
        if rhs == "1":
            return np.ones((len(test_individuals), 1), order="F")
        if data is None:
            msg = f"the formula is {rhs}, but no data is provided"
            raise ValueError(msg)
        if isinstance(data, str | Path):
            data = read_table(data)

        missing_individuals = set(test_individuals) - set(data.index)

        if missing_individuals:
            msg = f"The following individuals are missing in the provided data: {', '.join(missing_individuals)}"
            raise ValueError(msg)

        if not data.index.equals(pd.Index(test_individuals)):
            data = data.reindex(test_individuals)

        try:
            covariates = Formula(rhs).get_model_matrix(data, na_action="raise")
        except ValueError as e:
            error_msg = str(e)
            column_name = error_msg.split("`")[1]
            msg = f"`{column_name}` contains missing value. Which is unacceptable in the fixed effects."
            raise ValueError(msg) from None

        return np.asfortranarray(covariates)


class make_predictor:
    def __init__(self, train_bed: str | Path, chunk_size: int = 10000):
        self._train_bed = valid_path(train_bed, (".bed",))
        self._chunk_size = chunk_size

    def make(
        self,
        grms: dict[str, str | Path],
        params: str | Path | LinearMixedModel | LinearMixedModelParams,
    ):
        params = self._set_params(params)

        predictor = Predictor(self._train_bed, params)

        predictor._params = params
        predictor._random_effect_names = []
        predictor._lhs = params.lhs

        for name, grm in grms.items():
            self._set_grm(grm, predictor)
            predictor._random_effect_names.append(name)
        return predictor

    def _set_params(
        self,
        params: str | Path | LinearMixedModel | LinearMixedModelParams,
    ):
        if not isinstance(params, LinearMixedModelParams):
            return load_params(params)
        if isinstance(params, LinearMixedModelParams):
            return params

        msg = "Invalid parameter type. Please provide either a file path to parameters, a LinearMixedModel instance, or a LinearMixedModelParams object."
        raise ValueError(msg)

    def _set_grm(self, grm: str | Path, predictor: Predictor):
        if not hasattr(predictor, "_keep_alive"):
            predictor._keep_alive = []

        f = h5py.File(grm, "r")

        center = np.asfortranarray(f["center"][:])
        predictor.set_cross_grm(
            f.attrs["method"],
            center,
            f["scale_factor"][()],
            self._chunk_size,
        )

        f.close()
        predictor._keep_alive.extend([center])
