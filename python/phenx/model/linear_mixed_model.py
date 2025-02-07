import gc
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
        """
        Initialize the make_model class.

        Parameters
        ----------
        data : str | Path | pd.DataFrame
            Input data for the model. Can be a file path (str or Path) or a pandas DataFrame.
            If file path is provided, it will be read as a tab-separated file with the first
            column as index and NA values represented by ["NA", ".", "nan", "NaN"]. Other
            missing value is not support.
        categorical : list[str] | str | None, optional
            List of column names or single column name to be converted to categorical type.
            If None, no conversion is performed. Default is None.

        Raises
        ------
        ValueError
            If data is neither a file path nor a DataFrame.
        """
        if isinstance(data, (str | Path)):
            self.data = pd.read_csv(
                data, sep="\t", index_col=0, na_values=["NA", ".", "nan", "NaN"]
            )
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
        """
        Create a Linear Mixed Model from the specified formula and genetic relationship matrices.

        Parameters
        ----------
        formula : str
            A formula string specifying the model. The left-hand side specifies the response variable,
            and the right-hand side specifies the fixed effects. Example: "y ~ x1 + x2"
        grm : dict[str, pd.DataFrame | str | Path]
            A dictionary mapping random effect names to their corresponding genetic relationship matrices.
            The matrices can be provided as DataFrames or paths to files containing the matrices.

        Returns
        -------
        LinearMixedModel
            A LinearMixedModel object containing the response vector, design matrix, genetic relationship
            matrices, and random effect names.

        Raises
        ------
        ValueError
            If the fixed effects contain missing values.
        """
        formula = Formula(formula)
        data = self._clean_data(str(formula.lhs))
        data, self._grm, random_effect_names = self._load_grm(grm, data)

        try:
            model_matrix = formula.get_model_matrix(data, na_action="raise")
        except ValueError as e:
            msg = "fixed effects should not contain missing values"
            raise ValueError(msg) from e

        self._response = np.asfortranarray(model_matrix.lhs, dtype=np.float64)
        self._design_matrix = np.asfortranarray(model_matrix.rhs, dtype=np.float64)
        gc.collect()  # Clean up memory
        return LinearMixedModel(
            self._response,
            self._design_matrix,
            self._grm,
            random_effect_names,
        )

    def _clean_data(self, response_name: str) -> pd.DataFrame:
        """
        Remove rows containing missing values in the specified response column.

        Parameters
        ----------
        response_name : str
            Name of the response variable column in the DataFrame.

        Returns
        -------
        pd.DataFrame
            Dataframe with rows containing missing values in the response column removed.
        """
        return self.data.dropna(subset=[response_name])

    def _load_grm(
        self, grm: dict[str, pd.DataFrame | str | Path], data: pd.DataFrame
    ) -> tuple[np.ndarray, list[str]]:
        """
        Clean and align genetic relationship matrices (GRMs) with the data.

        Parameters
        ----------
        grm : dict[str, pd.DataFrame | str | Path]
            Dictionary mapping effect names to GRM matrices or paths to files.
        data : pd.DataFrame
            Dataframe containing the phenotypic and covariate data.

        Returns
        -------
        tuple[np.ndarray, list[str]]
            A tuple containing:
            - aligned_data: Dataframe with rows restricted to common indices
            - grm_cube: 3D numpy array of GRM matrices (samples x samples x effects)
            - random_effect_names: List of effect names

        Raises
        ------
        ValueError
            If no common indices are found between GRMs and data.
        """

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
