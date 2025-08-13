from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from gelexy._core import GRM, _BedReader


def make_grm(
    bed_file: str | Path,
    method: str = "add",
    chunk_size: int = 10000,
    save: bool = True,
) -> pd.DataFrame:
    """
    Construct a Genetic Relationship Matrix (GRM) from a BED genotype file.

    :param bed_file: Path to the BED genotype file.
    :param method: GRM calculation method, either 'add' (additive) or 'dom' (dominance). Default is 'add'.
    :param chunk_size: Number of SNPs to process per chunk. Default is 10,000.
    :param save: Whether to save the resulting GRM and metadata to an HDF5 file. Default is True.
    :raises ValueError: If an unsupported method is provided.
    :return: GRM as a pandas DataFrame with individuals as both index and columns.
    """
    bed_file = Path(bed_file)
    grm_maker = None
    method = method.lower()
    grm_maker = GRM(str(bed_file), chunk_size)
    add = True
    if method == "add":
        add = True
    elif method == "dom":
        add = False
    else:
        msg = f"Only support `add` and `dom` for method but got {method}."
        raise ValueError(msg)

    grm_arr = grm_maker.compute(add)
    grm_df = pd.DataFrame(
        grm_arr, index=grm_maker.individuals, columns=grm_maker.individuals
    )
    grm_df.index.name = "id"
    grm_df.columns.name = "id"

    if save:
        with h5py.File(f"{bed_file.with_suffix('')}.{method}.grm", "w") as f:
            f.create_dataset("grm", data=grm_arr)
            f.create_dataset("individuals", data=grm_maker.individuals)
            f.create_dataset("p_major", data=grm_maker.p_major)
            f.create_dataset("scale_factor", data=grm_maker.scale_factor)
            f.attrs["method"] = method
    return grm_df


def load_grm(
    path: str | Path,
    return_array: bool = False,
) -> pd.DataFrame | np.ndarray:
    """
    Load a Genetic Relationship Matrix (GRM) from an HDF5 file.

    :param path: Path to the HDF5 file containing the GRM.
    :param return_array: If True, return the GRM as a numpy ndarray; otherwise, return as a pandas DataFrame.
    :raises FileNotFoundError: If the specified GRM file does not exist.
    :return: The GRM as a pandas DataFrame or numpy ndarray, depending on `return_array`.
    """
    path = Path(path)
    if not path.exists():
        msg = f"GRM file {path} does not exist."
        raise FileNotFoundError(msg)

    grm = None
    with h5py.File(path, "r") as f:
        grm = f["grm"][:]
        individuals = f["individuals"].asstr()[:]
        grm = pd.DataFrame(grm, index=individuals, columns=individuals)
        grm.index.name = "id"
        grm.columns.name = "id"

    if return_array:
        return np.asfortranarray(grm)
    return grm


def load_genotype(bed_file: str, individuals: list[str] | None = None):
    if individuals is None:
        individuals = []
    reader = _BedReader(bed_file, int(1e10), individuals)
    genotype = pd.DataFrame(
        reader.read_chunk(),
        index=reader.individuals,
        columns=reader.snps,
    )
    genotype.index.name = "id"
    return genotype
