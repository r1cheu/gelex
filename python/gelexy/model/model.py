import gc

import numpy as np
import pandas as pd
from formulae import design_matrices
from scipy.sparse import csc_matrix

from gelexy import BayesAlphabet
from gelexy._core import Bayes, _BedReader, _sp_dense_dot

from .base import ModelMakerBase
from .formula_parser import FormulaParser


class make_bayes(ModelMakerBase):
    def make(
        self, formula: str, genotypes: dict[str], bayes_type: BayesAlphabet
    ):
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
        fparser = FormulaParser(formula, list(genotypes.keys()))
        data = self.data
        data = self._clear_data(data, fparser.response)

        for name, genotype in genotypes.items():
            genotypes[name] = load_genotype(genotype)

        data, design_matrix_genetic, dropped_ids, genotypes = (
            self._create_design_matrix_genetic(data, genotypes)
        )

        design_mat = design_matrices(fparser.common, data, na_action="error")

        model = Bayes(
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
                model.add_random_effect(
                    "(" + name + ")", np.asfortranarray(matrix.data)
                )

        for term in fparser.genetic_terms:
            genotypes = _sp_dense_dot(
                design_matrix_genetic, genotypes[term.genetic]
            )
            model.add_genetic_effect(
                term.name,
                genotypes,
                bayes_type,
            )

        model._dropped_ids = dropped_ids
        model._formula = formula
        gc.collect()
        return model

    def _create_design_matrix_genetic(
        self,
        data: pd.DataFrame,
        genotypes: dict[str, pd.DataFrame],
    ) -> tuple[np.ndarray, csc_matrix, list[str], dict[str, np.ndarray]]:
        data_ids = data["id"].to_numpy()
        grm_ids = next(iter(genotypes.values())).index.to_numpy()

        intersection = np.intersect1d(data_ids, grm_ids)
        dropped = np.setdiff1d(grm_ids, intersection)

        data_mask = np.isin(data_ids, intersection)
        data = data[data_mask]
        grm_mask = np.isin(grm_ids, intersection)

        for name, genotype in genotypes.items():
            genotypes[name] = genotype.loc[grm_mask, :]

        design_mat = make_design_matrix(
            data["id"].to_numpy(),
            next(iter(genotypes.values())).index.to_numpy(),
        )
        return data, design_mat, dropped, genotypes


def check_effect(formula: str, data: pd.DataFrame):
    try:
        model_matrix = design_matrices(formula, data, na_action="error")
    except ValueError as e:
        msg = "Common or Random effects columns contains missing values, which are unacceptable."
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


def load_genotype(bed_file: str):
    reader = _BedReader(bed_file, int(1e9))
    return pd.DataFrame(
        reader.read_chunk(),
        index=reader.individuals,
        columns=reader.snps,
    )


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
