from pathlib import Path

import h5py

from ..model import LinearMixedModel, LinearMixedModelParams


def load_params(model: str | Path | LinearMixedModel):
    params = None
    if isinstance(model, str | Path):
        model = Path(model)
        with h5py.File(model, "r") as f:
            beta = f["beta"][:]
            sigma = f["sigma"][:]
            dropped_individuals = f["dropped_individuals"].asstr()[:]
            params = LinearMixedModelParams(
                beta=beta, sigma=sigma, dropped_individuals=dropped_individuals
            )
            params._keep_alive = (beta, sigma, dropped_individuals)
            params.rhs = f.attrs["rhs"]

    if isinstance(model, LinearMixedModel):
        if not hasattr(model, "_dropped_individuals"):
            msg = "Model does not have dropped individuals. Are you get model from make_model.make?"
            raise RuntimeError(msg)

        if not hasattr(model, "_rhs"):
            msg = "Model does not have rhs. Are you get model from make_model.make?"
            raise RuntimeError(msg)

        params = LinearMixedModelParams(model, model._dropped_individuals)
        params.rhs = model._rhs

    return params
