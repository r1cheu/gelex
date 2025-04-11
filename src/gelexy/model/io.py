from pathlib import Path

import h5py
import numpy as np

from ..model import LinearMixedModel, LinearMixedModelParams


def load_params(model: str | Path | LinearMixedModel):
    params = None
    if isinstance(model, str | Path):
        model = Path(model)
        with h5py.File(model, "r") as f:
            params = LinearMixedModelParams(
                np.asfortranarray(f["beta"][:]),
                np.asfortranarray(f["sigma"][:]),
                np.asfortranarray(f["proj_y"][:]),
                f["dropped_individuals"].asstr()[:],
            )
            params.rhs = f.attrs["rhs"]
            params.lhs = f.attrs["lhs"]

    if isinstance(model, LinearMixedModel):
        try:
            params = LinearMixedModelParams(model, model._dropped_individuals)
            params.rhs = model._rhs
            params.lhs = model._lhs
        except AttributeError as e:
            msg = "Model does not have required attributes. Are you get model from make_model.make?"
            raise RuntimeError(msg) from e

    return params
