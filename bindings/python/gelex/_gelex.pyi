"""Python bindings for the gelex C++ library"""

import enum
from typing import Annotated

import numpy
from numpy.typing import NDArray


class GeneticMode(enum.Enum):
    A = 0

    D = 1

class GenotypeMethod(enum.Enum):
    StandardizeHWE = 0

    CenterHWE = 1

    Standardize = 2

    Center = 3

    OrthStandardizeHWE = 4

    OrthCenterHWE = 5

    OrthStandardize = 6

    OrthCenter = 7

    NOIAStandardize = 8

    NOIACenter = 9

def encode_inplace(genotypes: Annotated[NDArray[numpy.float64], dict(shape=(None, None))], mode: GeneticMode, method: GenotypeMethod) -> None: ...
