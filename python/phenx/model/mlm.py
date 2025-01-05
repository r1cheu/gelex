import logging
from datetime import datetime

import numpy as np
import pandas as pd
from formulaic import Formula
from numpy.typing import NDArray

from phenx import Phenotypes

from .._core import REMLLoop

logger = logging.getLogger(__name__)


class LinearMixedModel:
    def __init__(
        self,
        formula: str,
        rand: dict[str, pd.DataFrame],
        data: Phenotypes | pd.DataFrame,
        categorical: list[str] | str | None = None,
        use_double_precision: bool = True,
        drop_na: bool = True,
    ) -> None:
        """Initialize the Linear Mixed Model with optimized memory handling."""
        self._dtype = np.float64 if use_double_precision else np.float32
        self._setup_data(data, categorical, drop_na)
        self._process_formula(formula)
        self._setup_design_matrices()
        self._setup_random_effects(rand)

        # Initialize REML optimization
        self._reml = REMLLoop(
            self._response,
            self._design_matrix,
            self._z_index,
            self._rands,
            self._rand_names,
        )

        # Cache for expensive computations
        self._cache = {}

    def _setup_data(
        self,
        data: Phenotypes | pd.DataFrame,
        categorical: list[str] | str | None,
        drop_na: bool,
    ) -> None:
        if not isinstance(data, Phenotypes | pd.DataFrame):
            msg = "'data' must be a pandas DataFrame or `Phenotype`"
            raise TypeError(msg)

        data = data.data if isinstance(data, Phenotypes) else data
        self.data = with_categorical_cols(data, categorical)

        self._na_action = "drop" if drop_na else "raise"

    def _process_formula(self, formula: str) -> None:
        """Process formula with optimized matrix generation."""
        self._formula = formula
        formula_obj = Formula(formula)
        self._response_name = formula_obj.lhs

        # Generate model matrices efficiently
        self._model_matrix = formula_obj.get_model_matrix(
            self.data, na_action=self._na_action
        )
        self._full_model_matrix = formula_obj.get_model_matrix(
            self.data, na_action="ignore"
        )

    def _setup_design_matrices(self) -> None:
        """Set up design matrices with efficient memory layout."""
        if not hasattr(self._model_matrix, "lhs"):
            msg = "No outcome variable specified in formula"
            raise ValueError(msg)

        # Set up indices for efficient lookups
        self._index = self._model_matrix.lhs.index
        self._index_full = self._full_model_matrix.lhs.index
        self._index_nan = self._index_full.difference(self._index)

        # Convert matrices to efficient numpy arrays
        self._response = np.array(
            self._model_matrix.lhs.to_numpy(copy=True), dtype=self._dtype, order="F"
        )
        self._design_matrix = np.array(
            self._model_matrix.rhs.to_numpy(copy=True), dtype=self._dtype, order="F"
        )
        self._full_design = np.array(
            self._full_model_matrix.rhs.to_numpy(copy=True),
            dtype=self._dtype,
            order="F",
        )

        if np.isnan(self._full_design).any():
            msg = "Design matrix contains NaN values."
            raise ValueError(msg)

        self._log_setup_info()

    def _setup_random_effects(self, rands: dict[str, pd.DataFrame]) -> None:
        """Set up random effects with vectorized operations."""
        self._n_rand = len(rands)
        self._rand_names = list(rands)

        # Vectorized random effects processing
        try:
            self._rands = np.stack(
                [
                    mat.loc[self._index_full, self._index_full].to_numpy(copy=True)
                    for mat in rands.values()
                    if isinstance(mat, pd.DataFrame)
                ],
                axis=-1,
            ).astype(self._dtype, order="F")
        except KeyError as err:
            msg = "Phenotype and genotype samples do not match, possibly indicating more samples with only phenotype. Use `phenx.intersect` to get the intersection."
            raise KeyError(msg) from err

        if not self._rand_names:
            msg = "No random effects found. Ensure GRM is generated using px.grm and provided as dict of pandas DataFrames"
            raise ValueError(msg)

        self._formula += (
            f" + {' + '.join(f'u({name})' for name in self._rand_names)} + e"
        )
        self._z_index = self._index_full.get_indexer(self._index).astype(
            np.int64, order="F"
        )

    def _log_setup_info(self) -> None:
        """Log setup information efficiently."""
        logger.info(
            "Data summary:\nTotal samples: %d\nSamples with %s data: %d\nSamples to predict: %d",
            len(self._index_full),
            self._response_name,
            len(self._index),
            len(self._index_nan),
        )

    def fit(
        self,
        varcomp_prior: NDArray | None = None,
        method: str = "AI",
        em_init: bool = True,
        max_iter: int = 20,
        var_tolerance: float = 1e-8,
        verbose: bool = True,
    ) -> None:
        """Fit the model with optimized convergence."""
        if varcomp_prior is None:
            varcomp_prior = np.full(self._n_rand + 1, 1.0 / (self._n_rand + 1))

        logger.info(
            "Fitting model at %s\nFormula: %s",
            datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            self._formula,
        )

        self._reml.run(
            varcomp=varcomp_prior,
            method=method,
            em_init=em_init,
            max_iteration=max_iter,
            tolerance=var_tolerance,
            verbose=verbose,
        )

        # Clear cache after fitting
        self._cache.clear()

    @property
    def variance_component(self) -> NDArray:
        """Get variance components with caching using dict.get()."""
        return self._cache.get("variance_component") or self._cache.setdefault(
            "variance_component", self._reml.get_varcomp()
        )

    @property
    def beta(self) -> NDArray:
        """Get beta coefficients with caching using dict.get()."""
        return self._cache.get("beta") or self._cache.setdefault(
            "beta", self._reml.get_beta()
        )

    @property
    def blup(self) -> pd.DataFrame:
        """Get BLUPs with efficient DataFrame creation using dict.get()."""
        return self._cache.get("blup") or self._cache.setdefault(
            "blup",
            pd.DataFrame(
                data=self._reml.get_blup(),
                index=self._index_full,
                columns=self._rand_names,
            ),
        )

    @property
    def gebv(self) -> pd.DataFrame:
        """Get GEBVs with efficient computation using dict.get()."""
        return self._cache.get("gebv") or self._cache.setdefault(
            "gebv",
            pd.DataFrame(
                data=self._reml.get_gebv(self._full_design),
                index=self._index_full,
                columns=[self._response_name],
            ),
        )

    @property
    def gebv_pred(self) -> pd.DataFrame:
        """Get predicted GEBVs efficiently."""
        return self.gebv.loc[self._index_nan, :]

    @property
    def blup_pred(self) -> pd.DataFrame:
        """Get predicted BLUPs efficiently."""
        return self.blup.loc[self._index_nan, :]


def with_categorical_cols(data: pd.DataFrame, columns) -> pd.DataFrame:
    """Convert selected columns of a DataFrame to categorical type.

    It converts all object columns plus columns specified in the `columns` argument.
    """
    # Convert 'object' and explicitly asked columns to categorical.
    object_columns = list(data.select_dtypes("object").columns)
    to_convert = list(set(object_columns + listify(columns)))
    if to_convert:
        data = data.copy()  # don't modify original data frame
        data[to_convert] = data[to_convert].apply(lambda x: x.astype("category"))
    return data


def listify(obj):
    """Wrap all non-list or tuple objects in a list.

    Provides a simple way to accept flexible arguments.
    """
    if obj is None:
        return []

    return obj if isinstance(obj, (list | tuple | None)) else [obj]
