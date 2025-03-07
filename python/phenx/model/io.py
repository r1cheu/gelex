from pathlib import Path

import h5py
import numpy as np

from ..model import LinearMixedModel, LinearMixedModelParams


def load_params(model: str | Path | LinearMixedModel):
    params = None
    if isinstance(model, str | Path):
        model = Path(model)
        with h5py.File(model, "r") as f:
            beta = np.asfortranarray(f["beta"][:])
            sigma = np.asfortranarray(f["sigma"][:])
            proj_y = np.asfortranarray(f["proj_y"][:])
            dropped_individuals = f["dropped_individuals"].asstr()[:]

            params = LinearMixedModelParams(
                beta=beta,
                sigma=sigma,
                proj_y=proj_y,
                dropped_individuals=dropped_individuals,
            )

            params._keep_alive = (beta, sigma, proj_y)
            params.rhs = f.attrs["rhs"]
            params.lhs = f.attrs["lhs"]

    if isinstance(model, LinearMixedModel):
        if not hasattr(model, "_dropped_individuals"):
            msg = "Model does not have dropped individuals. Are you get model from make_model.make?"
            raise RuntimeError(msg)

        if not hasattr(model, "_rhs"):
            msg = "Model does not have rhs. Are you get model from make_model.make?"
            raise RuntimeError(msg)

        params = LinearMixedModelParams(model, model._dropped_individuals)
        params.rhs = model._rhs
        params.lhs = model._lhs

    return params
