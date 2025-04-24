from pathlib import Path

import pandas as pd


def read_table(path: str | Path):
    """
    Read a tab-separated file into a pandas DataFrame.

    Parameters
    ----------
    path : str or Path
        The file path to the tab-separated file. The first column is used as the index.
        Missing values are represented by "NA", ".", "nan", or "NaN".

    Returns
    -------
    pd.DataFrame
        A DataFrame with the first column as the index.
    """
    return pd.read_csv(path, sep="\t", na_values=["NA", ".", "nan", "NaN"])
