from pathlib import Path

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

    grm = grm_maker.compute()
    grm = pd.DataFrame(grm, index=grm_maker.individuals, columns=grm_maker.individuals)
    if save:
        grm.to_hdf(
            f"{bed_file.with_suffix('')}.{method}.grm",
            key="grm",
            mode="w",
        )
    return grm
