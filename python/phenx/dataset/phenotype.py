"""Phenotype class for read and process phenotype data."""

from pathlib import Path

import numpy as np
import pandas as pd

from .basetype import Basetypes


def list_to_str(input_list: list, num_cols: int) -> str:
    """
    Convert a list to a formatted string with a specified number of columns.

    Parameters
    ----------
    input_list : list
        Input list to be converted to string.
    num_cols : int
        Number of columns to display.

    Returns
    -------
    str
        Formatted string representation of the input list.
    """
    msg = ""
    max_length = max(len(str(item)) for item in input_list)
    for i, item in enumerate(input_list, start=1):
        msg += f"{item: <{max_length}}\t"
        if i % num_cols == 0:
            msg += "\n"
    return msg


class Phenotypes(Basetypes):
    """
    Phenotype class for reading and processing phenotype data.

    Attributes
    ----------
    num_phenotypes
    phenotypes
    """

    def _read_data(self, path: str | Path):
        """
        Read phenotype data from a file.

        Parameters
        ----------
        path : str | path
            Path to the phenotype data file.

        Returns
        -------
        pd.DataFrame
            Phenotype data as a pandas DataFrame.
        """
        return pd.read_csv(path, sep="\t", index_col=0)

    @property
    def num_phenotypes(self) -> int:
        """
        Return number of phenotypes in the dataset.

        Returns
        -------
        int
            Number of phenotypes in the dataset.
        """
        return self._data.shape[1]

    @property
    def phenotypes(self) -> list:
        """
        Return list of phenotype names.

        Returns
        -------
        list
            List of phenotype names.
        """
        return list(self._data.columns)

    def set_nan(self, input: float | int | list[int | float], seed=42) -> pd.DataFrame:
        """
        Set specified values in the DataFrame to NaN based on the input criteria.

        Parameters
        ----------
        input : float | int | list[int|float]
            The criteria for setting NaN values. If a float between 0 and 1, it represents
            the percentage of rows to set to NaN. If an int or float >= 1, it represents
            the number of rows to set to NaN. If a list, it contains the indices or labels
            of rows to set to NaN.
        seed : int, optional
            The seed for the random number generator, by default 42.

        Returns
        -------
        pd.DataFrame
            A copy of the DataFrame with specified values set to NaN.
        """

        rng = np.random.default_rng(seed)
        data_copy = self._data.copy()

        if isinstance(input, float) and 0 < input < 1:
            # Set a percentage of values to NaN
            num_nan = int(self.num_samples * input)
            nan_indices = rng.choice(self.num_samples, num_nan, replace=False)
            data_copy.iloc[nan_indices, :] = np.nan

        elif isinstance(input, int | float) and input >= 1:
            # Set a specific number of rows to NaN
            num_rows_to_nan = min(int(input), self.num_samples)
            nan_indices = rng.choice(self.num_samples, num_rows_to_nan, replace=False)
            data_copy.iloc[nan_indices, :] = np.nan

        elif isinstance(input, list):
            # Set values of the list (rows by index or labels) to NaN
            try:
                data_copy.loc[input, :] = np.nan  # Try to use labels
            except KeyError:
                data_copy.iloc[input, :] = np.nan  # Fallback to positional indexing

        return data_copy

    def __repr__(self):
        return (
            f"{self._path}:\n{self.num_samples} samples, {self.num_phenotypes} phenotypes\n"
            f"Samples:\n{list_to_str(self.samples[:5], 6)}...\n"
            f"Phenotypes:\n{list_to_str(self.phenotypes, 3)}"
        )

    def __getitem__(self, key):
        """
        Retrieve data from the DataFrame based on the provided key.

        Parameters
        ----------
        key : str or int
            The key to access the data. It can be a column label or row index.

        Returns
        -------
        pd.Series or pd.DataFrame
            The data corresponding to the provided key. If the key is a column label,
            a Series is returned. If the key is a row index, a DataFrame is returned.

        Raises
        ------
        KeyError
            If the key is not found in the DataFrame.
        """
        try:
            return self._data.loc[:, key]
        except KeyError:
            return self._data.loc[key, :]
