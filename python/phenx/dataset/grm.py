"""Genetic Relationship Matrix (GRM) computation."""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from phenx.utils.path import valid_path

from .._core import _Amat, _Amat_rbf, _Dmat, _Dmat_rbf
from .genotype import Genotypes


def grm(
    genotype: Genotypes | str | Path,
    encode: str = "add",
    kernel: str | None = None,
    save_grm: bool = False,
) -> pd.DataFrame:
    """
    Compute the Genetic Relationship Matrix (GRM) from genotype data.

    Parameters
    ----------
    genotype : Genotypes
        The genotype data to compute the GRM from. Can be a Genotypes object or a path to a `.h5` file.
    phenotype : pd.Series | str | Path | None, optional
        The phenotype data associated with the genotypes. Can be a pandas Series, a path to a file, or None.
    encode_method : str, optional
        The method to encode the genotype data. Default is "add".
    impute_method : str, optional
        The method to impute missing genotype data. Default is "median".
    save_grm : str | Path | None, optional
        The path to save the computed GRM. If None, the GRM is not saved.

    Returns
    -------
    pd.DataFrame
        The computed Genetic Relationship Matrix (GRM).

    Raises
    ------
    TypeError
        If the genotype is not a valid path to a `.h5` file or a Genotype object.
    """
    if isinstance(genotype, str | Path):
        genotype = valid_path(genotype, suffixes=(".h5",))
        logging.info("Loading genetic relationship matrix from %s", genotype)
        return pd.read_hdf(genotype, key="grm")

    if not isinstance(genotype, Genotypes):
        msg = "Genotype must be a path to a `.h5` file or a Genotype object."
        raise TypeError(msg)

    genotype_array = np.array(genotype, copy=True, dtype=np.float64)

    grm = None
    if encode == "add":
        if kernel is None:
            grm = _Amat(genotype_array)
        elif kernel == "rbf":
            bandwidth = np.sqrt(genotype_array.shape[0] / 2)
            grm = _Amat_rbf(genotype_array, bandwidth)
        else:
            msg = "Kernel method must be either None or 'rbf'."
            raise ValueError(msg)

    elif encode == "dom":
        if kernel is None:
            grm = _Dmat(genotype_array)
        elif kernel == "rbf":
            bandwidth = np.sqrt(genotype_array.shape[0] / 2)
            grm = _Dmat_rbf(genotype_array, bandwidth)
        else:
            msg = "Kernel method must be either None or 'rbf'."
            raise ValueError(msg)

    else:
        msg = "Genotype encoding method must be either 'add' or 'dom'."
        raise ValueError(msg)

    grm = pd.DataFrame(
        grm,
        index=genotype.data.columns,
        columns=genotype.data.columns,
    )

    if save_grm:
        grm.to_hdf(
            genotype.path.with_suffix(f".{encode}.h5"),
            key="grm",
        )

    return grm
