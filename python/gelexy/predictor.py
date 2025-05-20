# class Predictor(_Predictor):
#     def predict(
#         self,
#         test_bed: str | Path,
#         data: str | Path | pd.DataFrame | None = None,
#     ):
#         test_bed = valid_path(test_bed, (".bed",))
#         random_effects = self._compute_random_effects(test_bed)
#
#         prediction = random_effects.sum(
#             axis=1, keepdims=True
#         ) + self._compute_fixed_effects(
#             self._set_rhs(self.test_individuals, data)
#         )
#         return pd.DataFrame(
#             np.hstack([prediction, random_effects]),
#             index=self.test_individuals,
#             columns=[self._lhs, *self._random_effect_names],
#         )
#
#     def _set_rhs(
#         self,
#         test_individuals: list[str],
#         data: str | Path | pd.DataFrame | None = None,
#     ):
#         if self._rhs == "1":
#             return np.ones((len(test_individuals), 1), order="F")
#
#         if data is None:
#             msg = f"the formula is {self._rhs}, but no data is provided"
#             raise ValueError(msg)
#
#         if isinstance(data, str | Path):
#             data = read_table(data)
#
#         check_individuals(test_individuals, data)
#
#         data = data.loc[test_individuals]
#
#         common = check_effect(Formula(self._rhs), data).lhs
#
#         return np.asfortranarray(common, dtype=np.float64)
#
#
# def check_individuals(individuals: list[str], data: pd.DataFrame):
#     """Check if all specified individuals exist in the provided data.
#
#     Parameters
#     ----------
#     individuals : list[str]
#         List of individual identifiers to check
#     data : pd.DataFrame
#         DataFrame containing the individuals as index
#
#     Raises
#     ------
#     ValueError
#         If any individuals are missing from the data
#     """
#     missing_individuals = set(individuals) - set(data.index)
#     if missing_individuals:
#         msg = f"The following individuals are missing in the provided data: {', '.join(missing_individuals)}"
#         raise ValueError(msg)
#
#
# class make_predictor:
#     def __init__(self, train_bed: str | Path, chunk_size: int = 10000):
#         self._train_bed = valid_path(train_bed, (".bed",))
#         self._chunk_size = chunk_size
#
#     def make(
#         self,
#         grms: dict[str, str | Path],
#         params: str | Path | GBLUP | GBLUPParams,
#     ):
#         params = self._set_params(params)
#
#         predictor = Predictor(self._train_bed, params)
#
#         predictor._random_effect_names = []
#         predictor._lhs = params.lhs
#         predictor._rhs = params.rhs
#
#         for name, grm in grms.items():
#             self._set_cross_grm(grm, predictor)
#             predictor._random_effect_names.append(name)
#         return predictor
#
#     def _set_params(
#         self,
#         params: str | Path | GBLUP | GBLUPParams,
#     ):
#         if not isinstance(params, GBLUPParams):
#             return load_params(params)
#         if isinstance(params, GBLUPParams):
#             return params
#
#         msg = "Invalid parameter type. Please provide either a file path to parameters, a GBLUP instance, or a GBLUPParams object."
#         raise ValueError(msg)
#
#     def _set_cross_grm(self, grm: str | Path, predictor: Predictor):
#         with h5py.File(grm, "r") as f:
#             predictor.set_cross_grm(
#                 f.attrs["method"],
#                 np.asfortranarray(f["center"][:]),
#                 f["scale_factor"][()],
#                 self._chunk_size,
#             )
