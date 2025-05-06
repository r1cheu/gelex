import gc
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from formulae import design_matrices
from scipy.sparse import csc_matrix

from gelexy._core import _GBLUP
from gelexy.data.grm import load_grm

from ..data import read_table
from ..logger import setup_logger
from .formula_parser import FormulaParser


class GBLUP(_GBLUP):
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
            path = (
                f"{self._lhs}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.model"
            )

        with h5py.File(path, "w") as f:
            f.create_dataset("beta", data=self.beta)
            f.create_dataset("sigma", data=self.sigma)
            f.create_dataset("proj_y", data=self._proj_y)
            f.create_dataset("dropped_ids", data=self._dropped_ids)
            f.attrs["formula"] = self.formula()


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
            self.data = read_table(data)
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

    def make(self, formula: str, grms: dict) -> GBLUP:
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
        GBLUP
            A GBLUP object containing the response vector, design matrix, genetic relationship
            matrices, and random effect names.

        Raises
        ------
        ValueError
            If the fixed effects contain missing values.
        """
        fparser = FormulaParser(formula, list(grms.keys()))
        data = self.data
        data = self._clear_data(data, fparser.response)
        grms = self._load_grms(grms)

        data, design_matrix_genetic, dropped_ids, grms = (
            self._create_design_matrix_genetic(data, grms)
        )

        design_mat = design_matrices(fparser.common, data, na_action="error")

        model = GBLUP(
            fparser.format_common,
            np.asfortranarray(design_mat.response),
        )

        model.add_fixed_effect(
            list(design_mat.common.terms),
            self._get_fixed_levels(design_mat.common.terms),
            np.asfortranarray(design_mat.common.design_matrix),
        )

        if design_mat.group is not None:
            for name, matrix in design_mat.group.terms.items():
                model.add_group_effect(
                    "(" + name + ")", csc_matrix(matrix.data)
                )

        for term in fparser.genetic_terms:
            model.add_genetic_effect(
                term.name,
                design_matrix_genetic,
                np.asfortranarray(grms[term.genetic]),
            )

        for term in fparser.gxe_terms:
            dm = design_matrices("0+" + term.env, data, na_action="error")
            model.add_gxe_effect(
                term.name,
                design_matrix_genetic,
                np.asfortranarray(grms[term.genetic]),
                np.asfortranarray(dm.common),
            )

        model._dropped_ids = dropped_ids
        model._formula = formula
        gc.collect()
        return model

    def _clear_data(self, data: pd.DataFrame, response: str):
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

    @staticmethod
    def _load_grms(grms):
        for name, grm in grms.items():
            if isinstance(grm, str | Path):
                grms[name] = load_grm(grm)
            elif isinstance(grm, pd.DataFrame):
                grms[name] = grm
            else:
                msg = f"GRM for {name} must be a path to a file or a DataFrame."
                raise ValueError(msg)
        return grms

    def _create_design_matrix_genetic(
        self,
        data: pd.DataFrame,
        grms: dict[str, pd.DataFrame],
    ) -> tuple[np.ndarray, list[str]]:
        data_ids = data["id"].to_numpy()
        grm_ids = next(iter(grms.values())).index.to_numpy()

        intersection = np.intersect1d(data_ids, grm_ids)
        dropped = np.setdiff1d(grm_ids, intersection)

        data_mask = np.isin(data_ids, intersection)
        data = data[data_mask]
        grm_mask = np.isin(grm_ids, intersection)

        for name, grm in grms.items():
            grms[name] = grm.loc[grm_mask, grm_mask]

        design_mat = make_design_matrix(
            data["id"].to_numpy(), next(iter(grms.values())).index.to_numpy()
        )
        return data, design_mat, dropped, grms


def check_effect(formula: str, data: pd.DataFrame):
    try:
        model_matrix = design_matrices(formula, data, na_action="error")
    except ValueError as e:
        msg = "Common or Group effects columns contains missing values, which are unacceptable."
        raise ValueError(msg) from e
    return model_matrix


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


def make_design_matrix(obs_ids, grm_ids):
    id2idx = {id_: idx for idx, id_ in enumerate(grm_ids)}
    n_obs = len(obs_ids)

    rows = np.arange(n_obs)
    cols = np.empty(n_obs, dtype=np.int64)

    for i, id_ in enumerate(obs_ids):
        cols[i] = id2idx[id_]
    data = np.ones(n_obs, dtype=np.float64)
    return csc_matrix((data, (rows, cols)), shape=(n_obs, len(grm_ids)))
