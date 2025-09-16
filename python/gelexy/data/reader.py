"""
Data reading utilities for phenotype and plink fam files.

This module provides functions to read phenotype and .fam files commonly used in genetics.
It supports flexible handling of FID/IID columns and provides consistent index naming for downstream analysis.
"""

from pathlib import Path

import pandas as pd


def read_phenotype(
    path: str | Path, iid_only: bool = False, **kwargs
) -> pd.DataFrame:
    """
    Reads a phenotype file into a pandas DataFrame.

    :param path: Path to the phenotype file.
    :param iid_only: If True, expects only the 'IID' column and uses it as index.
                     If False, expects both 'FID' and 'IID' columns and creates a combined index.
    :param kwargs: Additional arguments passed to pandas.read_csv().
    :raises ValueError: If required columns are missing based on iid_only parameter.
    :return: DataFrame with individual IDs as index and phenotype data as columns.
    """
    path = Path(path).expanduser().resolve()

    phenotype = pd.read_csv(
        path,
        sep="\t",
        na_values=["NA", ".", "nan", "NaN"],
        dtype={"FID": str, "IID": str},
        **kwargs,
    )
    # we do not directly use the IID column as index to avoid issues with
    # duplicate IDs, will happen when multiple environment phenotypes are provided
    if iid_only:
        if "IID" not in phenotype.columns:
            msg = f"IID is needed, but not found in {path}. First three columns: {', '.join(phenotype.columns[:3])}"
            raise ValueError(msg)
        return phenotype

    if "FID" not in phenotype.columns or "IID" not in phenotype.columns:
        msg = f"FID and IID are needed, but not found in {path}. First three columns: {', '.join(phenotype.columns[:3])}"
        raise ValueError(msg)

    phenotype["FID_IID"] = phenotype["FID"] + "_" + phenotype["IID"]

    return phenotype


def read_fam(prefix: str, iid_only: bool = False) -> pd.Index:
    """
    Read a PLINK .fam file and return individual IDs as a pandas Index.

    :param prefix: prefix to the .fam file.
    :param iid_only: If True, return only the IID column as the index; otherwise, return FID+IID.
    :return: pandas.Index of individual IDs, named "IID" if iid_only else "FID_IID".
    """
    try:
        fam = pd.read_csv(
            str(prefix) + ".fam",
            sep=r"\s+",
            header=None,
            dtype={0: str, 1: str},
        )
    except FileNotFoundError as err:
        msg = f"{prefix}.fam not found. Are you sure the prefix is correct?"
        raise FileNotFoundError(msg) from err

    if iid_only:
        return pd.Index(fam.iloc[:, 1], name="IID")
    return pd.Index(fam.iloc[:, 0] + "_" + fam.iloc[:, 1], name="FID_IID")
