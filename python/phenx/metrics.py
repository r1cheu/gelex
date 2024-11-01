import logging
from abc import ABC, abstractmethod

import numpy as np
from scipy.stats import pearsonr


def tofloat(func):
    def wrapper(*args, **kwargs):  # numpydoc ignore=GL08
        return float(func(*args, **kwargs))

    return wrapper


class Metrics(ABC):
    def __init__(self):
        self._metric = []
        self._metric_name = self.__class__.__name__

    def update(self, true, pred):
        self._metric.append(self._compute(true, pred))

    @abstractmethod
    def _compute(self, true, pred):
        raise NotImplementedError

    def summary(self):
        logging.info(
            "%s\n    #run: %d, mean: %.4f, max: %.4f, min: %.4f",
            self._metric_name,
            len(self._metric),
            np.mean(self._metric),
            np.max(self._metric),
            np.min(self._metric),
        )

    def reset(self):
        self._metric = []


class Pearsonr(Metrics):
    @tofloat
    def _compute(self, true, pred):
        return pearsonr(true, pred)[0]


class RMSE(Metrics):
    @tofloat
    def _compute(self, true, pred):
        return np.sqrt(np.mean((true - pred) ** 2))


class MAE(Metrics):
    @tofloat
    def _compute(self, true, pred):
        return np.mean(np.abs(true - pred))


class R2(Metrics):
    @tofloat
    def _compute(self, true, pred):
        return 1 - np.sum((true - pred) ** 2) / np.sum((true - np.mean(true)) ** 2)


class MAPE(Metrics):
    @tofloat
    def _compute(self, true, pred):
        return np.mean(np.abs((true - pred) / true)) * 100.0

    def summary(self):
        logging.info(
            "%s\n    #run: %d, mean: %.4f%%, max: %.4f%%, min: %.4f%%",
            self._metric_name,
            len(self._metric),
            np.mean(self._metric),
            np.max(self._metric),
            np.min(self._metric),
        )
