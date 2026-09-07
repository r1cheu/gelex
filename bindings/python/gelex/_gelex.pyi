"""Python bindings for the gelex C++ library"""

from collections.abc import Sequence
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

class BinaryType(enum.Enum):
    float64 = 1

    float32 = 2

    int32 = 3

    uint8 = 4

class PayloadInfo:
    @property
    def identifier(self) -> str: ...

    @property
    def type(self) -> BinaryType: ...

    @property
    def shape(self) -> tuple[int, int]: ...

    def __repr__(self) -> str: ...

class PayloadWriterF64:
    """
    Handle to one reserved payload; append() adds one column (one draw) and write() stores the whole payload at once.
    """

    def append(self, column: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> None: ...

    def write(self, values: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='F', device='cpu', writable=False)]) -> None: ...

    @property
    def identifier(self) -> str: ...

    @property
    def rows(self) -> int: ...

class PayloadWriterF32:
    """
    Handle to one reserved payload; append() adds one column (one draw) and write() stores the whole payload at once.
    """

    def append(self, column: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> None: ...

    def write(self, values: Annotated[NDArray[numpy.float32], dict(shape=(None, None), order='F', device='cpu', writable=False)]) -> None: ...

    @property
    def identifier(self) -> str: ...

    @property
    def rows(self) -> int: ...

class PayloadWriterI32:
    """
    Handle to one reserved payload; append() adds one column (one draw) and write() stores the whole payload at once.
    """

    def append(self, column: Annotated[NDArray[numpy.int32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> None: ...

    def write(self, values: Annotated[NDArray[numpy.int32], dict(shape=(None, None), order='F', device='cpu', writable=False)]) -> None: ...

    @property
    def identifier(self) -> str: ...

    @property
    def rows(self) -> int: ...

class PayloadWriterU8:
    """
    Handle to one reserved payload; append() adds one column (one draw) and write() stores the whole payload at once.
    """

    def append(self, column: Annotated[NDArray[numpy.uint8], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> None: ...

    def write(self, values: Annotated[NDArray[numpy.uint8], dict(shape=(None, None), order='F', device='cpu', writable=False)]) -> None: ...

    @property
    def identifier(self) -> str: ...

    @property
    def rows(self) -> int: ...

class BinaryWriter:
    """
    Writer for gelex binary containers. Reserve payloads with a dtype and (rows, columns) shape, fill them column by column, then close() (or leave the with-block) to finalise the file.
    """

    def __init__(self, path: str) -> None: ...

    def reserve(self, identifier: str, type: BinaryType, shape: Sequence[int]) -> PayloadWriterF64 | PayloadWriterF32 | PayloadWriterI32 | PayloadWriterU8: ...

    def close(self) -> None: ...

    @property
    def is_open(self) -> bool: ...

    def __enter__(self) -> BinaryWriter: ...

    def __exit__(self, *args) -> None: ...

class BinaryReader:
    """
    Memory-mapped reader for gelex binary containers such as the MCMC .draws output. Payloads are exposed as read-only, column-major (rows, columns) NumPy views that alias the mapped file.
    """

    def __init__(self, path: str) -> None: ...

    def __len__(self) -> int: ...

    def __contains__(self, identifier: str) -> bool: ...

    def __getitem__(self, identifier: str) -> Annotated[NDArray, dict(writable=False)]: ...

    def info(self, identifier: str) -> PayloadInfo: ...

    def payloads(self) -> list[PayloadInfo]: ...

    def keys(self) -> list[str]: ...
