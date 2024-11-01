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
        self.data = self._read_data(self._path)

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
    def path(self) -> Path:
        """
        Return the path to the dataset.

        Returns
        -------
        Path
            The path to the dataset.
        """
        return self._path

    @property
    def head(self) -> pd.DataFrame:
        """
        Return the first few rows of the dataset.

        Returns
        -------
        pd.DataFrame
            The first few rows of the dataset.
        """
        return self.data.head()

    def __array__(self, dtype=None, copy=True) -> NDArray:
        """
        Convert the dataset to a NumPy array.

        Returns
        -------
        NDArray
            The dataset as a NumPy array.
        """
        return np.array(self.data, dtype=dtype, copy=copy)
