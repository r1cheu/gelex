# Copyright 2026 RuLei Chen
# SPDX-License-Identifier: Apache-2.0

"""Bridge from a gelex ``.draws`` container to ArviZ."""

from __future__ import annotations

import inspect
import os

import numpy as np

from ._gelex import CscReader, DenseReader, sparse_draws_path

__all__ = ["read_draws"]


def _marker_count(dense: DenseReader, sparse: CscReader | None) -> int | None:
    """Rows of the first genetic coefficients payload, or ``None``."""
    readers = [dense] if sparse is None else [sparse, dense]
    for reader in readers:
        for info in reader.payloads():
            parts = info.identifier.split("/")
            if len(parts) == 3 and parts[0] == "genetic" and parts[2] == "coefficients":
                return info.shape[0]
    return None


def read_draws(path: str | os.PathLike[str], *, include_markers: bool = False):
    """Read a ``.draws`` file (and its ``.csc`` companion) into ``arviz.InferenceData``.

    Every payload becomes a ``posterior`` variable with ``/`` replaced by
    ``.`` (xarray forbids ``/``) and a single chain: ``(1, draws)`` for
    scalar terms and ``(1, draws, rows)`` for vectors. Marker-level payloads
    (coefficients, assignments) are skipped unless ``include_markers`` is
    set, because their footprint is ``markers x draws``; the sparse ones are
    densified from the ``.csc`` companion when included.
    """
    import arviz as az

    dense_path = os.fspath(path)
    dense = DenseReader(dense_path)
    sparse_path = sparse_draws_path(dense_path)
    sparse = CscReader(sparse_path) if os.path.exists(sparse_path) else None

    markers = None if include_markers else _marker_count(dense, sparse)

    posterior: dict[str, np.ndarray] = {}
    dims: dict[str, list[str]] = {}

    def add(name: str, values: np.ndarray) -> None:
        name = name.replace("/", ".")  # xarray forbids "/" in names
        if values.shape[0] == 1:
            posterior[name] = np.array(values[0])[np.newaxis, :]
        else:
            posterior[name] = np.ascontiguousarray(values.T)[np.newaxis, :, :]
            dims[name] = [f"{name}_dim"]

    for info in dense.payloads():
        if info.shape[0] == markers:
            continue
        add(info.identifier, dense[info.identifier])

    if include_markers and sparse is not None:
        for info in sparse.payloads():
            add(info.identifier, sparse[info.identifier].toarray())

    if "posterior" in inspect.signature(az.from_dict).parameters:
        return az.from_dict(posterior=posterior, dims=dims)
    return az.from_dict({"posterior": posterior}, dims=dims)
