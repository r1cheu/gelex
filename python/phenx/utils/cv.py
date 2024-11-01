import logging

import numpy as np
import pandas as pd

from ..metrics import Metrics


class CrossValidation:
    def __init__(
        self, metrics: list[Metrics], n_splits: int = 5, seed: int = 42
    ) -> None:
        self._seed = seed
        self._rng = np.random.default_rng(seed)
        self._split = n_splits
        self._metrics = [metric() for metric in metrics]

    def test_index(self, num):
        indices = np.arange(num)
        self._rng.shuffle(indices)
        fold_sizes = np.full(self._split, num // self._split, dtype=int)
        fold_sizes[: num % self._split] += 1
        current = 0
        for fold_size in fold_sizes:
            start, stop = current, current + fold_size
            yield indices[start:stop]
            current = stop

    def evaluate(self, y_true: pd.DataFrame, y_pred: pd.DataFrame):
        logging.debug(
            "len of raw y_true %d, len of raw y_pred %d",
            len(y_true),
            len(y_pred),
        )
        y_true = y_true.dropna()
        y_pred = y_pred.loc[y_true.index, :]

        logging.debug(
            "len of y_true: %d, len of y_pred: %d, all match %s",
            len(y_true),
            len(y_pred),
            np.all(y_true.index == y_pred.index),
        )

        y_true = y_true.to_numpy().flatten()
        y_pred = y_pred.to_numpy().flatten()

        for metric in self._metrics:
            metric.update(y_true, y_pred)

    def summary(self):
        logging.info("Cross Validation Summary: Seed=%d", self._seed)
        for metric in self._metrics:
            metric.summary()

    def reset(self):
        for metric in self._metrics:
            metric.reset()
