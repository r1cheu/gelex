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

class DominanceCode(enum.Enum):
    Het = 0

    HWE = 1

    NOIA = 2

class Normalization(enum.Enum):
    None = 0

    Center = 1

    CenterScale = 2

class MomentBasis(enum.Enum):
    Empirical = 0

    Theoretical = 1

class LocusStats:
    @property
    def nA2A2(self) -> int: ...

    @property
    def nA1A2(self) -> int: ...

    @property
    def nA1A1(self) -> int: ...

    @property
    def n_missing(self) -> int: ...

    def n_nonmissing(self) -> int: ...

    def has_nonmissing(self) -> bool: ...

    def pA2A2(self) -> float: ...

    def pA1A2(self) -> float: ...

    def pA1A1(self) -> float: ...

    def A1freq(self) -> float: ...

class EncodingSpec:
    @property
    def effect(self) -> GeneticMode: ...

    @property
    def dominance_code(self) -> DominanceCode: ...

    @property
    def normalization(self) -> Normalization: ...

    @property
    def moment_basis(self) -> MomentBasis: ...

class LocusEncoding:
    @property
    def column_index(self) -> int: ...

    @property
    def marker_index(self) -> int: ...

    @property
    def stats(self) -> LocusStats: ...

    @property
    def code(self) -> list[float]: ...

    @property
    def mean(self) -> float: ...

    @property
    def var(self) -> float: ...

    @property
    def sd(self) -> float: ...

    @property
    def missing_encoded_value(self) -> float: ...

    @property
    def valid(self) -> bool: ...

class LociEncoding:
    @property
    def spec(self) -> EncodingSpec: ...

    @property
    def loci(self) -> list[LocusEncoding]: ...

def encode_inplace(genotypes: Annotated[NDArray[numpy.float64], dict(shape=(None, None))], effect: GeneticMode, method: GenotypeMethod, tol: float = 1e-12, marker_offset: int = 0) -> LociEncoding:
    """
    Encode a genotype matrix in place and return its LociEncoding.

    genotypes must be a writable float64 NumPy array in column-major
    (order='F') layout with shape (n_samples, n_markers); raw genotypes
    are 0/1/2 with NaN for missing. The array is modified in place (no
    copy); any other dtype/layout raises a TypeError.
    """
