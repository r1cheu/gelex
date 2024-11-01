"""Genotype class for read and process genotype data."""

from pathlib import Path

import numpy as np
import pandas as pd
from bed_reader import open_bed
from numpy.typing import NDArray

from phenx.utils import valid_path

from .basetype import Basetypes


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

    Attributes
    ----------
    ALLOWED_INPUT_SUFFIXES : tuple
        A tuple of allowed input file suffixes.
    """

    ALLOWED_INPUT_SUFFIXES = (".bed",)

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
            genotypes = np.asfortranarray(bed.read().T)
            np.nan_to_num(genotypes, copy=False, nan=1.0)
            snp_id = _concat_chrom_pos(bed.chromosome, bed.bp_position)
            samples = bed.iid

        return pd.DataFrame(genotypes, columns=samples, index=snp_id)

    def __repr__(self) -> str:
        mem_usage = self._get_memory_usage()
        num_loci = self.data.shape[0]
        num_samples = self.data.shape[1]
        file_path = self.path
        return (
            f"Genotypes(\n"
            f"    path='{file_path}',\n"
            f"    num_loci={num_loci},\n"
            f"    num_samples={num_samples},\n"
            f"    memory_usage={mem_usage:.2f} MB,\n"
            f")"
        )

    def _get_memory_usage(self) -> float:
        return self.data.memory_usage(deep=True).sum() / (1024**2)
