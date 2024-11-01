"""Genotype class for read and process genotype data."""

import logging
from abc import ABC, abstractmethod
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from phenx.utils import valid_path

from .._core import _hybrid, _hybrid_value  # noqa


class GenoProcessor(ABC):
    ALLOWED_METHODS = None

    def __init__(self, genotype, method):
        self._method = self._validate_method(method)
        self._index = genotype.data.index
        self._columns = genotype.data.columns
        self._dump_data = None

    def _check_columns(self, columns):
        if not all(columns == self._columns):
            msg = "loci in genotype and precompute data do not match."
            raise ValueError(msg)

    def _check_index(self, index):
        if not all(index == self._index):
            msg = "Samples in genotype and phenotype data do not match."
            raise ValueError(msg)

    def _validate_method(self, value: str) -> str:
        """
        Validate that the provided method is within the allowed methods.

        Parameters
        ----------
        value : str
            The method to validate.
        allowed_methods : list
            A list of allowed methods.

        Returns
        -------
        str
            The validated method.

        Raises
        ------
        ValueError
            If the method is not in the list of allowed methods.
        """
        if self.ALLOWED_METHODS is None:
            return value

        if value not in self.ALLOWED_METHODS:
            msg = f"`{value}` is not in {self.ALLOWED_METHODS}"
            raise ValueError(msg)

        return value

    @abstractmethod
    def run(self, genotype: NDArray):
        raise NotImplementedError


class GenoEncoder(GenoProcessor):
    ALLOWED_METHODS = ("add", "dom", "hybrid")

    def run(self, genotype: NDArray, phenotype: pd.Series | None):
        if self._method == "hybrid":
            self._hybrid_encode(genotype, phenotype)
        else:
            _encode(genotype, method=self._method)

    def _hybrid_encode(self, genotype: NDArray, phenotype: pd.Series) -> pd.DataFrame:
        """
        Perform hybrid encoding on the genotype data using the provided phenotype data.

        Parameters
        ----------
        genotype : NDArray
            The genotype data to be encoded.
        phenotype : pd.Series | Path | str | None
            The phenotype data used for hybrid encoding. Can be a pandas Series or a path to a file containing hybrid values.

        Returns
        -------
        pd.DataFrame
            The genotype data after hybrid encoding.

        Raises
        ------
        ValueError
            If phenotype data is not provided or if the samples in genotype and phenotype data do not match.
        """

        if phenotype is None:
            msg = "Phenotype data is required for hybrid encoding."
            raise ValueError(msg)

        self._check_index(phenotype.index)

        phenotype = np.array(phenotype, copy=True, dtype=np.float64)
        self._dump_data = _hybrid_value(genotype, phenotype)
        _hybrid(genotype, self._dump_data)

    def save(self, path: str | Path):
        path = Path(path)
        if self._method != "hybrid":
            return

        if self._dump_data is not None and not path.exists():
            pd.DataFrame(self._dump_data, columns=self._columns).to_hdf(
                path, key="encode"
            )
            logging.info("Encode values of each locus saved to %s", path)


class GenoFileEncoder(GenoProcessor):
    def run(self, genotype: NDArray, precompute: str | Path):
        self._dump_data = self._read_data(precompute)

        self._check_columns(self._dump_data.columns)

        values = np.array(self._dump_data, copy=True, dtype=np.float64)
        _hybrid(genotype, values)

    def _read_data(self, path: str | Path):
        path = valid_path(path, suffixes=(".h5",))
        return pd.read_hdf(path, key="encode")


class GenoImputer(GenoProcessor):
    ALLOWED_METHODS = ("mean", "median", "none")

    def run(self, genotype: NDArray):
        if self._method == "none":
            return
        self._dump_data = _impute(genotype, method=self._method)

    def save(self, path: str | Path):
        path = Path(path)
        if self._method == "none":
            return
        if self._dump_data is None:
            msg = "No data to save. invoke `run` first."
            raise ValueError(msg)
        if not path.exists():
            pd.DataFrame(self._dump_data, index=self._columns).to_hdf(
                path, key="impute"
            )
            logging.info("Imputed values of each locus saved to %s", path)


class GenoFileImputer(GenoProcessor):
    def run(self, genotype: NDArray, precompute: str | Path):
        self._dump_data = self._read_data(precompute)

        self._check_columns(self._dump_data.columns)
        _value_impute(genotype, self._dump_data)

    def _read_data(self, path: str | Path):
        path = valid_path(path, suffixes=(".h5",))
        return pd.read_hdf(path, key="impute")
