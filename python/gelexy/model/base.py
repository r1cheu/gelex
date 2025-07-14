from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix

from .._logger import setup_logger
from ..data import read_pheno


class ModelMakerBase:
    def __init__(
        self,
        data: str | Path | pd.DataFrame,
        categorical: list[str] | str | None = None,
    ):
        """
        Initialize the base model with phenotype data.

        :param data: Path to a file, pandas DataFrame, or pandas Series containing phenotype data.
        :param categorical: List of column names or a single column name to treat as categorical, or None.
        :raises ValueError: If data is not a valid type or does not contain an 'id' column.
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

        self._logger = setup_logger(__name__)

    def _dropna(self, data: pd.DataFrame, response: str) -> pd.DataFrame:
        """
        Drops rows from the DataFrame where the response column contains missing values.

        :param data: Input DataFrame to process.
        :param response: Name of the response column to check for missing values.
        :return: DataFrame with rows containing missing values in the response column removed.
        """
        if data[response].hasnans:
            self._logger.info(
                "Missing values detected in `%s`. These entries will be dropped.",
                response,
            )
            return data.dropna(subset=[response])
        return data


def get_fixed_levels(fixed_effect: dict):
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


def make_design_matrix(
    phenotype_ids: list[str], genotype_ids: list[str]
) -> csc_matrix:
    """
    Constructs a sparse design matrix mapping phenotype IDs to genotype IDs.

    :param phenotype_ids: List of phenotype identifiers (observations).
    :param genotype_ids: List of genotype identifiers (features/columns).
    :return: A csc_matrix of shape (len(phenotype_ids), len(genotype_ids)), where each row has a 1 at the column corresponding to the genotype ID for that phenotype.
    """
    id2idx = {id_: idx for idx, id_ in enumerate(genotype_ids)}
    n_obs = len(phenotype_ids)

    rows = np.arange(n_obs)
    cols = np.empty(n_obs, dtype=np.int64)

    for i, id_ in enumerate(phenotype_ids):
        cols[i] = id2idx[id_]
    data = np.ones(n_obs, dtype=np.float64)
    return csc_matrix((data, (rows, cols)), shape=(n_obs, len(genotype_ids)))
