import numpy as np
from formulae.matrices import DesignMatrices

from gelexy import BayesModel, load_genotype
from gelexy._core import MCMCResult, _BayesPredictor


class BayesPredictor(_BayesPredictor):
    def __init__(self, model: BayesModel, result: MCMCResult):
        super().__init__(model, result)
        self._design_matrix: DesignMatrices = model._design_matrix

    def predict(self, data, genotypes: dict):
        fixed = np.array(
            self._design_matrix.common.evaluate_new_data(data), order="F"
        )
        random = np.empty(shape=(1, 1, 1), order="F")
        if self._design_matrix.group is not None:
            random = np.array(
                self._design_matrix.group.evaluate_new_data(data),
                order="F",
            )
        if fixed.ndim == 1:
            fixed = fixed[:, np.newaxis]
        if random.ndim == 2:
            random = random[:, :, np.newaxis]

        genotype = np.array(
            np.stack(
                [
                    load_genotype(genotype)
                    for key, genotype in genotypes.items()
                ],
                axis=-1,
            ),
            order="F",
        )

        return self._predict(fixed, random, genotype)
