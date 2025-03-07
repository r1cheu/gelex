import gc
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from formulaic import Formula

from .._chenx import _LinearMixedModel, _LinearMixedModelParams
from ..data import load_grm, read_table
from ..logger import setup_logger


class LinearMixedModel(_LinearMixedModel):
    @property
    def U(self):
        return pd.DataFrame(
            self._U, index=self._individuals, columns=self.random_effect_names
        )

    def save(self, path: str | Path | None = None):
        """
        Save the model to a file.

        Parameters
        ----------
        path : str | Path | None
            Path to the file where the model will be saved.
        """
        if path is None:
            path = f"{self._lhs}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.model"

        with h5py.File(path, "w") as f:
            f.create_dataset("beta", data=self.beta)
            f.create_dataset("sigma", data=self.sigma)
            f.create_dataset("proj_y", data=self._proj_y)
            f.create_dataset("dropped_individuals", data=self._dropped_individuals)
            f.attrs["rhs"] = self._rhs
            f.attrs["lhs"] = self._lhs


class LinearMixedModelParams(_LinearMixedModelParams):
    pass  # inherit from _LinearMixedModelParams, for the ability to extend class in python


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
            Input data for the model. See phenx.read_table for details on accepted formats.
        categorical : list[str] | str | None, optional
            List of column names or single column name to be converted to categorical type.
            If None, no conversion is performed. Default is None.

        Raises
        ------
        ValueError
            If data is neither a file path nor a DataFrame.
        """
        if isinstance(data, (str | Path)):
            self.data = read_table(data)
        elif isinstance(data, pd.DataFrame):
            self.data = data
        elif isinstance(data, pd.Series):
            self.data = pd.DataFrame(data)
        else:
            msg = "data must be a path to a file, DataFrame or Series object."
            raise ValueError(msg)

        if categorical:
            self.data = with_categorical_cols(self.data, categorical)
        self._logger = setup_logger(__name__)

    def make(self, formula: str, grm: dict[str, str | Path]) -> LinearMixedModel:
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
        self._lhs = str(formula.lhs)
        rhs = str(formula.rhs)

        data = self.data

        data, dropped_individuals = self._clear_data(data)

        grm_cube, data, random_effect_names, dropped_individuals = self._load_grm(
            grm, data, dropped_individuals
        )

        try:
            model_matrix = formula.get_model_matrix(data, na_action="raise")
        except ValueError as e:
            error_msg = str(e)
            column_name = error_msg.split("`")[1]
            msg = f"`{column_name}` contains missing value. Which is unacceptable in the fixed effects."
            raise ValueError(msg) from None

        response = np.asfortranarray(model_matrix.lhs, dtype=np.float64)
        design_matrix = np.asfortranarray(model_matrix.rhs, dtype=np.float64)

        model = LinearMixedModel(
            response,
            design_matrix,
            grm_cube,
            random_effect_names,
        )
        model._keep_alive = (
            response,
            design_matrix,
            grm_cube,
        )  # need keep the memory address valid

        model._dropped_individuals = list(dropped_individuals)
        model._rhs = rhs
        model._lhs = self._lhs
        gc.collect()  # grm is quite large, we need to collect the garbage

        return model

    def _clear_data(self, data: pd.DataFrame):
        if data[self._lhs].hasnans:
            dropped_individuals = data.index[data[self._lhs].isna()]
            self._logger.info(
                "Missing values detected in `%s`. These entries will be dropped.",
                self._lhs,
            )
            return data.dropna(subset=[self._lhs]), dropped_individuals
        return data, pd.Index([])

    def _load_grm(
        self,
        grm: dict[str, pd.DataFrame | str | Path],
        data: pd.DataFrame,
        drop_individuals: set[str],
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

        data_index = data.index
        grm_matrices = []
        random_effect_names = []

        for effect_name, grm_value in grm.items():
            matrix = (
                load_grm(grm_value)
                if isinstance(grm_value, (str | Path))
                else grm_value
            )

            if not data_index.isin(matrix.index).all():
                msg = f"Some individuals in the `{self._lhs}` are not present in the GRM. Are you sure you are using the correct GRM?"
                raise ValueError(msg)

            if len(matrix.index) != len(data_index):
                drop_individuals = drop_individuals.union(
                    matrix.index.difference(data_index)
                )
                matrix = matrix.drop(index=drop_individuals, columns=drop_individuals)

            grm_matrices.append(matrix)
            random_effect_names.append(str(effect_name))

        grm_cube = np.asfortranarray(np.stack(grm_matrices, axis=-1))

        return (
            grm_cube,
            data.loc[matrix.index],
            random_effect_names,
            drop_individuals,
        )


def with_categorical_cols(data: pd.DataFrame, columns) -> pd.DataFrame:
    """Convert selected columns of a DataFrame to categorical type.

    It converts all object columns plus columns specified in the `columns` argument.
    """
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
