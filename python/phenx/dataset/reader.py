"""Contains help functions for reading and writing hdf."""

import pandas as pd


def save_hdf(data: pd.DataFrame, path: str, key: str, meta: dict | None = None):
    """
    Save a pandas DataFrame to an HDF5 file.

    Parameters
    ----------
    data : pd.DataFrame
        The DataFrame to be saved.
    path : str
        The file path where the HDF5 file will be saved.
    key : str
        The key under which the DataFrame will be stored in the HDF5 file.
    meta : dict | None, optional
        A dictionary of metadata to be stored with the DataFrame. Default is None.
    """
    with pd.HDFStore(path) as store:
        store.put(key, data)
        if meta is not None:
            for k, v in meta.items():
                store.get_storer(key).attrs[k] = v


def read_hdf(path: str, key: str, meta: bool = True) -> pd.DataFrame:
    """
    Read a pandas DataFrame from an HDF5 file.

    Parameters
    ----------
    path : str
        The file path from which the HDF5 file will be read.
    key : str
        The key under which the DataFrame is stored in the HDF5 file.
    meta : bool, optional
        If True, return the metadata associated with the DataFrame. Default is True.

    Returns
    -------
    pd.DataFrame
        The DataFrame read from the HDF5 file.
    dict, optional
        A dictionary of metadata associated with the DataFrame, if meta is True.
    """
    with pd.HDFStore(path) as store:
        local_metadata = dict(store.get_storer(key).attrs)
        if meta:
            return store[key], local_metadata
    return store[key]
