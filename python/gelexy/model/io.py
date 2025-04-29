from pathlib import Path

import h5py
import numpy as np

from ..model import GBLUP, GBLUPParams


def load_params(model: str | Path | GBLUP):
    params = None
    if isinstance(model, str | Path):
        model = Path(model)
        with h5py.File(model, "r") as f:
            params = GBLUPParams(
                np.asfortranarray(f["beta"][:]),
                np.asfortranarray(f["sigma"][:]),
                np.asfortranarray(f["proj_y"][:]),
                f["dropped_ids"].asstr()[:],
            )
            params.rhs = f.attrs["rhs"]
            params.lhs = f.attrs["lhs"]

    if isinstance(model, GBLUP):
        try:
            params = GBLUPParams(model, model._dropped_ids)
            params.rhs = model._rhs
            params.lhs = model._lhs
        except AttributeError as e:
            msg = "Model does not have required attributes. Are you get model from make_model.make?"
            raise RuntimeError(msg) from e

    return params
