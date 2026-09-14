"""Python bindings for the gelex C++ library"""

from collections.abc import Mapping, Sequence
import enum
from typing import Annotated

import numpy
from numpy.typing import NDArray
import scipy.sparse


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

class BinaryType(enum.Enum):
    float64 = 1

    float32 = 2

    uint8 = 4

class MatrixHeader:
    @property
    def identifier(self) -> str: ...

    @property
    def type(self) -> BinaryType: ...

    @property
    def shape(self) -> tuple[int, int]: ...

    def __repr__(self) -> str: ...

class DenseStreamF64:
    """
    Handle to one reserved matrix; append() adds one column (one draw) and write() stores the whole matrix at once.
    """

    def append(self, column: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> None: ...

    def write(self, values: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='F', device='cpu', writable=False)]) -> None: ...

    @property
    def identifier(self) -> str: ...

class DenseStreamF32:
    """
    Handle to one reserved matrix; append() adds one column (one draw) and write() stores the whole matrix at once.
    """

    def append(self, column: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> None: ...

    def write(self, values: Annotated[NDArray[numpy.float32], dict(shape=(None, None), order='F', device='cpu', writable=False)]) -> None: ...

    @property
    def identifier(self) -> str: ...

class DenseStreamU8:
    """
    Handle to one reserved matrix; append() adds one column (one draw) and write() stores the whole matrix at once.
    """

    def append(self, column: Annotated[NDArray[numpy.uint8], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> None: ...

    def write(self, values: Annotated[NDArray[numpy.uint8], dict(shape=(None, None), order='F', device='cpu', writable=False)]) -> None: ...

    @property
    def identifier(self) -> str: ...

class DenseWriter:
    """
    Writer for gelex dense containers. Reserve matrices with a dtype and (rows, columns) shape, fill them column by column, then close() (or leave the with-block) to publish the file; every matrix must be complete. An unclosed writer discards its output.
    """

    def __init__(self, path: str) -> None: ...

    def reserve(self, identifier: str, type: BinaryType, shape: Sequence[int]) -> DenseStreamF64 | DenseStreamF32 | DenseStreamU8: ...

    def close(self) -> None: ...

    @property
    def is_open(self) -> bool: ...

    def __enter__(self) -> DenseWriter: ...

    def __exit__(self, *args) -> None: ...

class CscStreamF64:
    """
    Handle to one reserved sparse matrix; append() adds one dense column (one draw) and stores its non-zero entries.
    """

    def append(self, column: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> None: ...

class CscStreamF32:
    """
    Handle to one reserved sparse matrix; append() adds one dense column (one draw) and stores its non-zero entries.
    """

    def append(self, column: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> None: ...

class CscStreamU8:
    """
    Handle to one reserved sparse matrix; append() adds one dense column (one draw) and stores its non-zero entries.
    """

    def append(self, column: Annotated[NDArray[numpy.uint8], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> None: ...

class CscWriter:
    """
    Writer for gelex CSC containers such as the MCMC .draws.csc output. Reserve matrices with a dtype and (rows, columns) shape, append dense columns whose zeros are dropped, then close() (or leave the with-block) to publish the file; every matrix must be complete. An unclosed writer discards its output.
    """

    def __init__(self, path: str) -> None: ...

    def reserve(self, identifier: str, type: BinaryType, shape: Sequence[int]) -> CscStreamF64 | CscStreamF32 | CscStreamU8: ...

    def close(self) -> None: ...

    @property
    def is_open(self) -> bool: ...

    def __enter__(self) -> CscWriter: ...

    def __exit__(self, *args) -> None: ...

class DenseReader:
    """
    Memory-mapped reader for gelex binary containers such as the MCMC .draws output. Payloads are exposed as read-only, column-major (rows, columns) NumPy views that alias the mapped file.
    """

    def __init__(self, path: str) -> None: ...

    def __len__(self) -> int: ...

    def __contains__(self, identifier: str) -> bool: ...

    def __getitem__(self, identifier: str) -> Annotated[NDArray, dict(writable=False)]: ...

    def info(self, identifier: str) -> MatrixHeader: ...

    def payloads(self) -> list[MatrixHeader]: ...

    def keys(self) -> list[str]: ...

class CscReader:
    """
    Memory-mapped reader for gelex CSC containers such as the MCMC .draws.csc output. Matrices are returned as scipy.sparse.csc_matrix copies of shape (rows, columns).
    """

    def __init__(self, path: str) -> None: ...

    def __len__(self) -> int: ...

    def __contains__(self, identifier: str) -> bool: ...

    def __getitem__(self, identifier: str) -> scipy.sparse.csc_matrix[float] | scipy.sparse.csc_matrix[float] | scipy.sparse.csc_matrix[int]: ...

    def info(self, identifier: str) -> MatrixHeader: ...

    def nnz(self, identifier: str) -> int: ...

    def payloads(self) -> list[MatrixHeader]: ...

    def keys(self) -> list[str]: ...

def sparse_draws_path(draws_path: str) -> str:
    """Path of the CSC companion written next to a .draws file."""

class Bed:
    """
    A PLINK1 .bed dataset with its .fam/.bim metadata. gather() narrows and reorders the samples; every design built afterwards follows that order.
    """

    @property
    def num_samples(self) -> int: ...

    @property
    def num_markers(self) -> int: ...

    @property
    def sample_ids(self) -> list[tuple[str, str]]:
        """(FID, IID) of every sample in the current order."""

    @property
    def marker_ids(self) -> list[str]: ...

    def gather(self, samples: Sequence[tuple[str, str]]) -> None:
        """
        Restrict the dataset to the given (FID, IID) pairs, in that order. Every pair must exist in the .fam file.
        """

def open_bed(bfile: str) -> Bed:
    """Open `bfile`.bed/.bim/.fam."""

class GeneticProjection:
    """
    One genetic mode's encoded view of a design's genotypes: each marker is a 4-entry lookup table over the raw BED codes.
    """

    @property
    def num_samples(self) -> int: ...

    @property
    def num_markers(self) -> int: ...

    @property
    def snp_luts(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='F')]:
        """
        (4, markers) lookup tables indexed by raw BED code (A1A1, missing, A1A2, A2A2).
        """

def gebv(projection: GeneticProjection, coefficients: Annotated[NDArray[numpy.float64], dict(shape=(None,), writable=False)]) -> Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]:
    """
    Genomic values of every sample under `projection` for one vector of marker coefficients (one posterior draw).
    """

class GeneticDesign:
    """
    The genotypes of one Bed's current samples together with a projection per genetic mode.
    """

    @property
    def num_samples(self) -> int: ...

    @property
    def num_markers(self) -> int: ...

    @property
    def modes(self) -> list[GeneticMode]: ...

    @property
    def a1_frequency(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')]: ...

    def contains(self, mode: GeneticMode) -> bool: ...

    def projection(self, mode: GeneticMode) -> GeneticProjection: ...

def make_genetic_design(bed: Bed, luts: Mapping[GeneticMode, Annotated[NDArray[numpy.float64], dict(shape=(4, None), order='F')]]) -> GeneticDesign:
    """
    Build the design from the lookup tables of a trained model, as returned by load_snp_luts(). The Bed must carry the same markers in the same order and allele orientation as the training data.
    """

def load_snp_luts(path: str) -> dict[GeneticMode, Annotated[NDArray[numpy.float64], dict(shape=(4, None), order='F')]]:
    """Read a .snplut file into {GeneticMode: (4, markers) array}."""

def write_snp_luts(path: str, luts: Mapping[GeneticMode, Annotated[NDArray[numpy.float64], dict(shape=(4, None), order='F')]]) -> None: ...
