from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from ._chenx import _Predictor
from .data import load_grm
from .model import LinearMixedModel, LinearMixedModelParams
from .model.io import load_params


class Predictor(_Predictor):
    def predict(self, test_bed: str):
        u = self._predict(test_bed)
        return pd.DataFrame(u, index=self.test_individuals)


class make_predictor:
    def __init__(self, train_bed: str, chunk_size: int):
        self._train_bed = train_bed
        self._chunk_size = chunk_size

    def make(
        self,
        grms: dict[str, str | Path],
        params: str | Path | LinearMixedModel | LinearMixedModelParams,
    ):
        params = self._set_params(params)
        dropped_individuals = params.dropped_individuals.copy()
        predictor = Predictor(self._train_bed, params)
        for _, grm in grms.items():
            self._set_grm(grm, dropped_individuals, predictor)
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

    def _set_grm(
        self, grm: str | Path, dropped_individuals: list[str], predictor: Predictor
    ):
        if not hasattr(predictor, "_keep_alive"):
            predictor._keep_alive = []

        f = h5py.File(grm, "r")
        grm = load_grm(grm, return_array=True, dropped_individual=dropped_individuals)
        predictor.set_grm(grm)

        center = np.asfortranarray(f["center"][:])
        predictor.set_cross_grm(
            f.attrs["method"],
            center,
            f["scale_factor"][()],
            self._chunk_size,
        )

        f.close()

        predictor._keep_alive.extend([grm, center])
