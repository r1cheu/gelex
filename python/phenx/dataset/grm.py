"""Genetic Relationship Matrix (GRM) computation."""

import logging
from pathlib import Path

import pandas as pd

from phenx.utils import valid_path

from .._core import _grm
from .genotype import Genotypes
from .reader import save_hdf

logging.basicConfig(format="%(message)s", level=logging.INFO)


def grm(
    genotype: Genotypes,
    phenotype: pd.Series | str | Path | None = None,
    encode_method: str = "add",
    impute_method: str = "mean",
    save_grm: str | Path | None = None,
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
        The method to impute missing genotype data. Default is "mean".
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

    if isinstance(genotype, Genotypes):
        genotype.encode = encode_method
        genotype.impute = impute_method

    elif isinstance(genotype, str | Path):
        genotype = valid_path(genotype, suffixes=(".h5"))
        logging.info("Loading genetic relationship matrix from %s", genotype)
        return pd.read_hdf(genotype, key="grm")
    else:
        msg = "Genotype must be a path to a `.h5` file or a Genotype object."
        raise TypeError(msg)

    if encode_method == "hybrid":
        grm = pd.DataFrame(
            _grm(genotype(phenotype), encode_method=encode_method),
            index=genotype.samples,
            columns=genotype.samples,
        )
    else:
        grm = pd.DataFrame(
            _grm(genotype(), encode_method=encode_method),
            index=genotype.samples,
            columns=genotype.samples,
        )

    if save_grm is not None:
        save_grm = Path(save_grm)

        if save_grm.suffix != ".h5":
            warning = "The file extension is not `.h5`. The GRM will be saved with `.h5` extension."
            logging.warning(warning)
            save_grm = save_grm.with_suffix(f".{encode_method}.h5")

        save_hdf(save_grm, grm, key="grm")
        logging.info("GRM saved to %s", save_grm)

    return grm
