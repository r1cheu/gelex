import h5py

from .linear_mixed_model import LinearMixedModel


def save_model(path: str, model: LinearMixedModel):
    with h5py.File(path, "w") as f:
        f.create_group("model")
        group = f["model"]
        group.create_dataset("beta", model.beta)
        group.create_dataset("sigma", model.sigma)
        group.create_dataset("individual_names", model.individual_names)
