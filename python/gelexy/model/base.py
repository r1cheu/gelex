from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix

from ..data import read_pheno
from ..logger import setup_logger


class ModelMakerBase:
    def __init__(
        self,
        data: str | Path | pd.DataFrame,
        categorical: list[str] | str | None = None,
    ):
        """
        Initialize the model maker base class.

        Parameters
        ----------
        data : str | Path | pd.DataFrame
            Input data for the model. See gelexy.read_table for details on accepted formats.
        categorical : list[str] | str | None, optional
            List of column names or single column name to be converted to categorical type.
            If None, no conversion is performed. Default is None.

        Raises
        ------
        ValueError
            If data is neither a file path nor a DataFrame.
        """
        if isinstance(data, (str | Path)):
            self.data = read_pheno(data)
        elif isinstance(data, pd.DataFrame):
            self.data = data
        elif isinstance(data, pd.Series):
            self.data = pd.DataFrame(data)
        else:
            msg = "data must be a path to a file, DataFrame or Series object."
            raise ValueError(msg)
        if "id" not in self.data.columns:
            msg = "data must contain an 'id' column."
            raise ValueError(msg)

        if categorical:
            self.data = with_categorical_cols(self.data, categorical)

        self.aligner = None
        self._logger = setup_logger(__name__)

    def _dropna(self, data: pd.DataFrame, response: str):
        if data[response].hasnans:
            self._logger.info(
                "Missing values detected in `%s`. These entries will be dropped.",
                response,
            )
            return data.dropna(subset=[response])
        return data

    @staticmethod
    def _get_fixed_levels(fixed_effect: dict):
        levels = []
        for name, term in fixed_effect.items():
            if hasattr(term, "levels"):
                levels.extend([f"{name}_{level}" for level in term.levels])
            else:
                levels.append(name)
        return levels


def with_categorical_cols(data: pd.DataFrame, columns) -> pd.DataFrame:
    """Convert selected columns of a DataFrame to categorical type.

    It converts all object columns plus columns specified in the `columns` argument.
    """
    object_columns = list(data.select_dtypes("object").columns)
    to_convert = list(set(object_columns + listify(columns)))
    if to_convert:
        data[to_convert] = data[to_convert].apply(
            lambda x: x.astype("category")
        )
    return data


def listify(obj):
    """Wrap all non-list or tuple objects in a list.

    Provides a simple way to accept flexible arguments.
    """
    if obj is None:
        return []

    return obj if isinstance(obj, (list | tuple | None)) else [obj]


def make_design_matrix(phenotype_ids, genotype_ids):
    id2idx = {id_: idx for idx, id_ in enumerate(genotype_ids)}
    n_obs = len(phenotype_ids)

    rows = np.arange(n_obs)
    cols = np.empty(n_obs, dtype=np.int64)

    for i, id_ in enumerate(phenotype_ids):
        cols[i] = id2idx[id_]
    data = np.ones(n_obs, dtype=np.float64)
    return csc_matrix((data, (rows, cols)), shape=(n_obs, len(genotype_ids)))


def create_design_matrix_genetic(
    data: pd.DataFrame,
    genotypes: dict[pd.DataFrame],
    symmetric: bool = True,
) -> tuple[np.ndarray, list[str]]:
    data_ids = data["id"].to_numpy()
    genotypes_ids = next(iter(genotypes.values())).index.to_numpy()

    intersection = np.intersect1d(data_ids, genotypes_ids)
    dropped = np.setdiff1d(genotypes_ids, intersection)

    data_mask = np.isin(data_ids, intersection)
    data = data[data_mask]
    genotype_mask = np.isin(genotypes_ids, intersection)

    for name, genotype in genotypes.items():
        if symmetric:
            genotypes[name] = genotype.loc[genotype_mask, genotype_mask]
        else:
            genotypes[name] = genotype.loc[genotype_mask, :]

    design_mat = make_design_matrix(
        data["id"].to_numpy(),
        next(iter(genotypes.values())).index.to_numpy(),
    )
    return data, design_mat, dropped, genotypes
