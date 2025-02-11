from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from .._chenx import add_grm, dom_grm


def make_grm(
    bed_file: str | Path,
    method: str = "add",
    chunk_size: int = 10000,
    save: bool = True,
) -> pd.DataFrame:
    """
    Compute the Genetic Relationship Matrix (GRM) from a BED file.

    Parameters
    ----------
    bed_file : str | Path
        Path to the BED(PLINK) file containing SNP data.
    method : str, optional
        Method to use for GRM computation. Must be either 'add' (additive) or 'dom' (dominance).
        Default is 'add'.
    chunk_size : int | bool, optional
        Number of SNPs to process in each chunk. If False, the entire dataset is processed at once.
        Default is 10000.
    save : bool, optional
        Save the computed GRM to a file. Default is True.

    Returns
    -------
    np.ndarray
        The computed Genetic Relationship Matrix (GRM).

    Raises
    ------
    ValueError
        If the method is not 'add' or 'dom'.
    """
    bed_file = Path(bed_file)
    grm_maker = None
    method = method.lower()
    if method == "add":
        grm_maker = add_grm(str(bed_file), chunk_size)
    elif method == "dom":
        grm_maker = dom_grm(str(bed_file), chunk_size)
    else:
        msg = f"Only support `add` and `dom` for method but got {method}."
        raise ValueError(msg)

    grm_arr = grm_maker.compute()
    grm_df = pd.DataFrame(
        grm_arr, index=grm_maker.individuals, columns=grm_maker.individuals
    )
    if save:
        with h5py.File(f"{bed_file.with_suffix('')}.{method}.grm", "w") as f:
            f.create_dataset("grm", data=grm_arr)
            f.create_dataset("individuals", data=grm_maker.individuals)
            f.create_dataset("center", data=grm_maker.center)
            f.create_dataset("scale_factor", data=grm_maker.scale_factor)
            f.attrs["method"] = method
    return grm_df


def load_grm(
    grm_path: str | Path, return_array: bool = False
) -> pd.DataFrame | np.ndarray:
    """
    Load GRM from HDF5 file.

    Parameters
    ----------
    grm_path : str | Path
        Path to GRM HDF5 file.
    return_array : bool, optional
        Return as NumPy array if True, else pandas DataFrame. Default False.

    Returns
    -------
    pd.DataFrame | np.ndarray
        Loaded GRM data.

    Raises
    ------
    FileNotFoundError
        If GRM file not found.
    """
    grm_path = Path(grm_path)
    if not grm_path.exists():
        msg = f"GRM file {grm_path} does not exist."
        raise FileNotFoundError(msg)

    with h5py.File(grm_path, "r") as f:
        grm = np.asfortranarray(f["grm"][:])
        if not return_array:
            individuals = f["individuals"][:]
            return pd.DataFrame(grm, index=individuals, columns=individuals)

    return grm
