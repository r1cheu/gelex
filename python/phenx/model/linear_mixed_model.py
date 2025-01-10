from pathlib import Path

import numpy as np
import pandas as pd
from formulaic import Formula

from .._chenx import LinearMixedModel


class make_model:
    def __init__(
        self,
        data: str | Path | pd.DataFrame,
        categorical: list[str] | str | None = None,
    ):
        if isinstance(data, (str | Path)):
            self.data = pd.read_csv(data, sep="\t", index_col=0)
        elif isinstance(data, pd.DataFrame):
            self.data = data
        else:
            msg = "data must be a path to a file or a DataFrame"
            raise ValueError(msg)
        if categorical:
            self.data = with_categorical_cols(self.data, categorical)

    def make(
        self, formula: str, grm: dict[str, pd.DataFrame | str | Path]
    ) -> LinearMixedModel:
        formula = Formula(formula)
        data = self._clean_data(str(formula.lhs))
        data, grm, random_effect_names = self._clean_grm(grm, data)

        try:
            model_matrix = formula.get_model_matrix(data, na_action="raise")
        except ValueError as e:
            msg = "fixed effects should not contain missing values"
            raise ValueError(msg) from e

        response = np.asfortranarray(model_matrix.lhs, dtype=np.float64)
        design_matrix = np.asfortranarray(model_matrix.rhs, dtype=np.float64)

        return LinearMixedModel(
            response,
            design_matrix,
            grm,
            random_effect_names,
        )

    def _clean_data(self, response_name: str) -> pd.DataFrame:
        return self.data.dropna(subset=[response_name])

    def _clean_grm(
        self, grm: dict[str, pd.DataFrame | str | Path], data: pd.DataFrame
    ) -> tuple[np.ndarray, list[str]]:
        grm_matrices = []
        random_effect_names = []
        common_index = None

        for effect_name, grm_value in grm.items():
            # Load GRM matrix if path is provided
            matrix = (
                pd.read_hdf(grm_value)
                if isinstance(grm_value, (str | Path))
                else grm_value
            )

            # Find common indices between GRM and data
            if common_index is None:
                common_index = matrix.index.intersection(data.index)
            # Align matrix with common indices
            aligned_matrix = matrix.loc[common_index, common_index]
            grm_matrices.append(aligned_matrix)
            random_effect_names.append(str(effect_name))

        aligned_data = data.loc[common_index]
        grm_cube = np.stack(grm_matrices, axis=-1)

        if not grm_cube.flags["F_CONTIGUOUS"]:
            grm_cube = np.asfortranarray(grm_cube)

        return aligned_data, grm_cube, random_effect_names


def with_categorical_cols(data: pd.DataFrame, columns) -> pd.DataFrame:
    """Convert selected columns of a DataFrame to categorical type.

    It converts all object columns plus columns specified in the `columns` argument.
    """
    # Convert 'object' and explicitly asked columns to categorical.
    object_columns = list(data.select_dtypes("object").columns)
    to_convert = list(set(object_columns + listify(columns)))
    if to_convert:
        data[to_convert] = data[to_convert].apply(lambda x: x.astype("category"))
    return data


def listify(obj):
    """Wrap all non-list or tuple objects in a list.

    Provides a simple way to accept flexible arguments.
    """
    if obj is None:
        return []

    return obj if isinstance(obj, (list | tuple | None)) else [obj]
