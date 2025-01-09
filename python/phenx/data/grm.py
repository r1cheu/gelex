from pathlib import Path

import numpy as np
import pandas as pd
from bed_reader import open_bed

from .._chenx import add_grm, add_grm_chunk, dom_grm, dom_grm_chunk


# TODO: current we depend on the bed_reader package to read the bed file.
# Thus, the chunk function is less efficient than it could be. When the
# pgnelib is decoupled from plink2, we can utilize it and write c++ code
# to read the bed file and compute the GRM in chunks.
# invoke rust from c++ is kind of ugly.
def make_grm(
    bed_file: str | Path,
    method: str = "add",
    chunk_size: int | bool = 10000,
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
    bed = open_bed(bed_file)
    grm = (
        np.zeros((bed.shape[0], bed.shape[0]), order="F", dtype=np.float64)
        if chunk_size
        else None
    )

    chunk_func_map = {"add": (add_grm_chunk, add_grm), "dom": (dom_grm_chunk, dom_grm)}

    if method not in chunk_func_map:
        msg = "method must be either 'add' or 'dom'."
        raise ValueError(msg)

    chunk_func, oneshot_func = chunk_func_map[method]

    if chunk_size:
        grm = _chunk_grm(chunk_func, bed, grm, chunk_size)
    else:
        grm = oneshot_func(bed.read(dtype=np.float64))

    grm = pd.DataFrame(grm, index=bed.iid, columns=bed.iid)
    if save:
        grm.to_hdf(
            f"{bed_file.with_suffix('')}.{method}.grm",
            key="grm",
            mode="w",
            format="table",
        )
    return grm


def _load_grm(
    grm_file: str | Path,
    samples_col: list[str] | None = None,
    samples_row: list[str] | None = None,
) -> np.ndarray:
    """
    Load a Genetic Relationship Matrix (GRM) from a file.
    Usually users don't need to call this function directly.

    Parameters
    ----------
    grm_file : str | Path
        Path to the file containing the GRM.
    samples_col : list[str], optional
        List of column sample IDs to load from the GRM. If None, loads all columns.
    samples_row : list[str], optional
        List of row sample IDs to load from the GRM. If None, loads all rows.

    Returns
    -------
    np.ndarray
        The loaded Genetic Relationship Matrix (GRM) as a NumPy array.
    """
    grm_file = Path(grm_file)

    return pd.read_hdf(
        grm_file,
        columns=samples_col,
        where=("index in %r" % samples_row) if samples_row else None,  # noqa: UP031
    ).to_numpy()


def _chunk_grm(
    chunk_func: callable, bed: open_bed, grm: np.ndarray, chunk_size: int
) -> np.ndarray:
    """
    Compute GRM in chunks.

    Parameters
    ----------
    chunk_func : callable
        Function to compute GRM for a chunk of SNP data
    bed : open_bed
        BED file object containing SNP data
    grm : np.ndarray
        GRM matrix to be updated
    chunk_size : int
        Number of SNPs to process in each chunk

    Returns
    -------
    np.ndarray
        Updated GRM matrix after processing all chunks
    """
    num_snps = bed.shape[1]
    for start in range(0, num_snps, chunk_size):  # split the SNP data into chunks
        end = min(start + chunk_size, num_snps)
        snp_data = bed.read(np.s_[:, start:end], dtype=np.float64)
        chunk_func(snp_data, grm)
    return grm
