# Copyright 2026 RuLei Chen
# SPDX-License-Identifier: Apache-2.0

"""Bridge from a gelex ``.draws`` container to ArviZ."""

from __future__ import annotations

import inspect
import os

import numpy as np

from ._gelex import DenseReader

__all__ = ["read_draws"]


def read_draws(path: str | os.PathLike[str], *, include_markers: bool = False):
    """Read a ``.draws`` file into an ``arviz.InferenceData``.

    Every payload becomes a ``posterior`` variable with ``/`` replaced by
    ``.`` (xarray forbids ``/``) and a single chain: ``(1, draws)`` for
    scalar terms and ``(1, draws, rows)`` for vectors. Marker-sized payloads
    (coefficients, assignments, unpooled marker variances) are skipped unless
    ``include_markers`` is set, because their footprint is ``markers x draws``.
    """
    import arviz as az

    reader = DenseReader(os.fspath(path))
    payloads = reader.payloads()

    markers = None
    if not include_markers:
        for info in payloads:
            parts = info.identifier.split("/")
            if len(parts) == 3 and parts[0] == "genetic" and parts[2] == "coefficients":
                markers = info.shape[0]
                break

    posterior: dict[str, np.ndarray] = {}
    dims: dict[str, list[str]] = {}
    for info in payloads:
        if info.shape[0] == markers:
            continue
        payload = reader[info.identifier]
        name = info.identifier.replace("/", ".")  # xarray forbids "/" in names
        if payload.shape[0] == 1:
            posterior[name] = np.array(payload[0])[np.newaxis, :]
        else:
            posterior[name] = np.ascontiguousarray(payload.T)[np.newaxis, :, :]
            dims[name] = [f"{name}_dim"]

    if "posterior" in inspect.signature(az.from_dict).parameters:
        return az.from_dict(posterior=posterior, dims=dims)
    return az.from_dict({"posterior": posterior}, dims=dims)
