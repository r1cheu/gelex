"""Genotype class for read and process genotype data."""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from bed_reader import open_bed
from numpy.typing import NDArray

from phenx.dataset.reader import save_hdf
from phenx.utils import valid_path

from .._core import _encode, _hybrid, _hybrid_value, _impute
from .basetype import Basetypes

logging.basicConfig(format="%(message)s", level=logging.INFO)


class Genotypes(Basetypes):
    """
    A class to handle genotype datasets in the phenx project.

    The Genotypes class extends the Basetypes class and provides additional
    functionalities specific to genotype data, including encoding and imputing
    missing values.

    Parameters
    ----------
    path : str, Path
        The path to the genotype data file.
    encode : str, optional
        The encoding method to be used (default is "add").
    impute : str, optional
        The imputation method to be used (default is "mean").

    Attributes
    ----------
    ALLOWED_ENCODE_METHODS : tuple
        A tuple of allowed encoding methods.
    ALLOWED_IMPUTE_METHODS : tuple
        A tuple of allowed imputation methods.
    ALLOWED_INPUT_SUFFIXES : tuple
        A tuple of allowed input file suffixes.
    encode
    impute
    """

    ALLOWED_ENCODE_METHODS = ("add", "dom", "hybrid")
    ALLOWED_IMPUTE_METHODS = ("mean", "median")
    ALLOWED_INPUT_SUFFIXES = (".bed",)

    def __init__(self, path: str | Path, encode: str = "add", impute: str = "mean"):
        super().__init__(path)
        self._encode = self._validate_method(
            encode, self.ALLOWED_ENCODE_METHODS, "encode_method"
        )
        self._impute = self._validate_method(
            impute, self.ALLOWED_IMPUTE_METHODS, "impute_method"
        )

    def _read_data(self, path: str | Path) -> pd.DataFrame:
        """
        Read plink bed file and return genotype data as pd.DataFrame.

        Parameters
        ----------
        path : str | Path
            The path to the genotype data file.

        Returns
        -------
        pd.DataFrame
            The genotype data. Rows are samples and columns are SNPs loci.
        """

        def _concat_chrom_pos(chrom: NDArray, pos: NDArray, sep: str = "_") -> NDArray:
            """
            Concatenate chromosome and position arrays into a single array with a separator.

            Parameters
            ----------
            chrom : NDArray
                An array of chromosome values.
            pos : NDArray
                An array of position values.
            sep : str, optional
                A separator string to use between chromosome and position. Defaults to "_".

            Returns
            -------
            NDArray
                An array with concatenated chromosome and position values.
            """
            return np.char.add(np.char.add(chrom.astype(str), sep), pos.astype(str))

        path = valid_path(path, suffixes=(".bed",))

        with open_bed(path) as bed:
            genotypes = np.asfortranarray(bed.read())
            snp_id = _concat_chrom_pos(bed.chromosome, bed.bp_position)
            samples = bed.iid

        return pd.DataFrame(genotypes, columns=snp_id, index=samples)

    @property
    def encode(self):
        """
        Get the current encoding method.

        Returns
        -------
        str
            The current encoding method.
        """
        return self._encode

    @encode.setter
    def encode(self, value: str):
        """
        Set the encoding method.

        Parameters
        ----------
        value : str
            The encoding method to be set.

        Raises
        ------
        ValueError
            If the encoding method is not in the allowed encoding methods.
        """
        self._encode = self._validate_method(
            value, self.ALLOWED_ENCODE_METHODS, "encode_method"
        )

    @property
    def impute(self):
        """
        Get the current imputation method.

        Returns
        -------
        str
            The current imputation method.
        """
        return self._impute

    @impute.setter
    def impute(self, value: str):
        """
        Set the imputation method.

        Parameters
        ----------
        value : str
            The imputation method to be set.

        Raises
        ------
        ValueError
            If the imputation method is not in the allowed imputation methods.
        """
        self._impute = self._validate_method(
            value, self.ALLOWED_IMPUTE_METHODS, "impute_method"
        )

    def __repr__(self):
        return (
            f"Genotype(path={self._path}, encode={self._encode}, impute={self._impute})"
        )

    def __call__(self, phenotype: pd.Series | None = None) -> pd.DataFrame:
        """
        Process the genotype data using the specified encoding and imputation methods.

        Parameters
        ----------
        phenotype : pd.Series, optional
            An optional pandas Series representing the phenotype data.

        Returns
        -------
        pd.DataFrame
            A DataFrame containing the processed genotype data.
        """
        logging.info(
            "Processing genotype using `%s` encode and impute nan using `%s`",
            self._encode,
            self._impute,
        )
        genotype = np.array(self._data.to_numpy(), copy=True, dtype=np.float64)
        if self._encode == "hybrid":
            self._hybrid_encode(genotype, phenotype)
        else:
            _encode(genotype, method=self._encode)

        _impute(genotype, method=self._impute)

        return genotype

    def _validate_method(
        self, value: str, allowed_methods: tuple[str], method_name: str
    ) -> str:
        """
        Validate that the provided method is within the allowed methods.

        Parameters
        ----------
        value : str
            The method to validate.
        allowed_methods : list
            A list of allowed methods.
        method_name : str
            The name of the method being validated.

        Returns
        -------
        str
            The validated method.

        Raises
        ------
        ValueError
            If the method is not in the list of allowed methods.
        """
        if value not in allowed_methods:
            msg = f"{method_name} must be one of {allowed_methods}"
            raise ValueError(msg)
        return value

    def _hybrid_encode(
        self, genotype: NDArray, phenotype: pd.Series | Path | str | None
    ) -> pd.DataFrame:
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
            msg = "Phenotype data is required for hybrid encoding. Either provide a path to a file contain hybird_values(.h5) or a pandas Series contain Phenotype data(see doc of Phenotypes)."
            raise ValueError(msg)

        hybrid_values = None

        if isinstance(phenotype, str | Path):
            phenotype = valid_path(phenotype, suffixes=(".h5",))
            hybrid_values = pd.read_hdf(phenotype, key="hybird_values")

        elif isinstance(phenotype, pd.Series):
            msg = "samples in genotype and phenotype data do not match. Try use phenotype = phenotype & genotype to take intersect of both"
            if len(phenotype) != self._data.shape[0]:
                raise ValueError(msg)

            if not all(phenotype.index == self._data.index):
                raise ValueError(msg)
            phenotype = np.array(phenotype.to_numpy(), copy=True, dtype=np.float64)
            hybrid_values = _hybrid_value(genotype, phenotype)
            save_path = self._path.with_suffix(".hybird_values.h5")
            save_hdf(
                pd.DataFrame(hybrid_values, columns=self._data.columns),
                save_path,
                key="hybird_values",
            )
            logging.info("Hybrid values saved to %s", save_path)

        _hybrid(genotype, hybrid_values)
        return genotype
