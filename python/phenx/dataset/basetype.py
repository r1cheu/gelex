"""Base Abstract class for handling genotype and phenotype in the phenx."""

from abc import ABC, abstractmethod
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from phenx.utils import valid_path


class Basetypes(ABC):
    """
    Abstract base class for handling datasets in the phenx project.

    The Basetypes class provides a template for dataset operations, including methods for reading data, accessing data properties, and converting data to different formats. Subclasses must implement the _read_data method to define how data is read from a given path.

    Parameters
    ----------
    path : str | path
        The path to the data file.

    Attributes
    ----------
    data
    head
    samples
    num_samples
    """

    ALLOWED_INPUT_SUFFIXES = None

    def __init__(self, path: str | Path):
        """
        Initialize the Basetypes instance.

        Parameters
        ----------
        path : str | Path
            The path to the data file.
        """
        self._path = valid_path(path, suffixes=self.ALLOWED_INPUT_SUFFIXES)
        self._data = self._read_data(self._path)

    @abstractmethod
    def _read_data(self, path) -> pd.DataFrame:
        """
        Read data from the given path.

        Parameters
        ----------
        path : str | Path
            The path to the data file.

        Returns
        -------
        pd.DataFrame
            The data read from the file.
        """
        raise NotImplementedError

    @property
    def data(self) -> pd.DataFrame:
        """
        Return a copy of the dataset.

        Returns
        -------
        pd.DataFrame
            A copy of the dataset.
        """
        return self._data.copy()

    @property
    def head(self) -> pd.DataFrame:
        """
        Return the first few rows of the dataset.

        Returns
        -------
        pd.DataFrame
            The first few rows of the dataset.
        """
        return self._data.head()

    @property
    def samples(self) -> list:
        """
        Return a list of sample indices from the dataset.

        Returns
        -------
        list
            A list of sample indices.
        """
        return list(self._data.index)

    @property
    def num_samples(self) -> int:
        """
        Return the number of samples in the dataset.

        Returns
        -------
        int
            The number of samples in the dataset.
        """
        return len(self.samples)

    def to_array(self) -> NDArray:
        """
        Convert the dataset to a NumPy array.

        Returns
        -------
        NDArray
            The dataset as a NumPy array.
        """
        return np.array(self._data)

    def __array__(self) -> NDArray:
        """
        Convert the dataset to a NumPy array (used for interoperability with NumPy).

        Returns
        -------
        NDArray
            The dataset as a NumPy array.
        """
        return self.to_array()

    def __and__(self, other) -> tuple:
        """
        Perform an intersection of samples with another `Basetypes`.

        Both `Basetypes` instances are modified in place and return.

        Parameters
        ----------
        other : Basetypes
            Another instance of Basetypes to intersect with.

        Returns
        -------
        tuple
            A tuple containing the two Basetypes instances with intersected data.
        """
        if not isinstance(other, Basetypes):
            msg = f"unsupported operand type(s) for &: 'Basetypes' and '{type(other)}'"
            raise TypeError(msg)
        intersection = self._data.index.intersection(other._data.index)
        intersection = sorted(intersection)
        self._data = self._data.loc[intersection]
        other._data = other._data.loc[intersection]
        return self, other
