# Copyright 2026 RuLei Chen
# SPDX-License-Identifier: Apache-2.0

import gc
from pathlib import Path

import gelexy
import numpy as np
import pytest

# 5 samples x 4 markers as additive dosages; NaN is missing and marker 2 is
# monomorphic, so a trained model leaves its lookup table all-zero.
DOSAGES = np.array(
    [
        [0.0, 2.0, 1.0, 0.0],
        [1.0, 2.0, 1.0, np.nan],
        [2.0, 1.0, 1.0, 2.0],
        [0.0, 0.0, 1.0, 1.0],
        [1.0, 2.0, 1.0, 2.0],
    ]
)

# Raw BED code per dosage: 0 -> A2A2 (3), 1 -> A1A2 (2), 2 -> A1A1 (0), NaN -> 1.
CODE_OF = {0.0: 3, 1.0: 2, 2.0: 0}


def write_plink(prefix: Path, dosages: np.ndarray) -> None:
    samples, markers = dosages.shape
    with open(f"{prefix}.fam", "w") as fam:
        for i in range(samples):
            fam.write(f"fam{i} id{i} 0 0 1 -9\n")
    with open(f"{prefix}.bim", "w") as bim:
        for j in range(markers):
            bim.write(f"1 snp{j} 0 {j + 1} A G\n")
    payload = bytearray([0x6C, 0x1B, 0x01])
    for j in range(markers):
        column = bytearray((samples + 3) // 4)
        for i in range(samples):
            value = dosages[i, j]
            code = 1 if np.isnan(value) else CODE_OF[value]
            column[i // 4] |= code << (2 * (i % 4))
        payload.extend(column)
    Path(f"{prefix}.bed").write_bytes(payload)


def centred_luts(dosages: np.ndarray) -> np.ndarray:
    """(4, markers) centred additive lookup tables indexed by raw BED code."""
    mean = np.nanmean(dosages, axis=0)
    luts = np.stack([2.0 - mean, np.zeros_like(mean), 1.0 - mean, -mean])
    luts[:, np.nanstd(dosages, axis=0) == 0] = 0.0
    return np.asfortranarray(luts)


def centred_design(dosages: np.ndarray) -> np.ndarray:
    """Dense centred additive design with missing values imputed to the mean."""
    mean = np.nanmean(dosages, axis=0)
    design = np.where(np.isnan(dosages), mean, dosages) - mean
    design[:, np.nanstd(dosages, axis=0) == 0] = 0.0
    return design


@pytest.fixture
def bfile(tmp_path: Path) -> str:
    prefix = tmp_path / "toy"
    write_plink(prefix, DOSAGES)
    return str(prefix)


@pytest.fixture
def snplut(tmp_path: Path) -> str:
    path = str(tmp_path / "toy.snplut")
    gelexy.write_snp_luts(path, {gelexy.GeneticMode.A: centred_luts(DOSAGES)})
    return path


def test_bed_exposes_axes_and_gathers(bfile: str):
    bed = gelexy.open_bed(bfile)
    assert bed.num_samples == 5
    assert bed.num_markers == 4
    assert bed.marker_ids == [f"snp{j}" for j in range(4)]
    assert bed.sample_ids == [(f"fam{i}", f"id{i}") for i in range(5)]

    picked = [("fam4", "id4"), ("fam1", "id1"), ("fam0", "id0")]
    bed.gather(picked)
    assert bed.num_samples == 3
    assert bed.sample_ids == picked

    with pytest.raises(RuntimeError):
        bed.gather([("fam9", "id9")])
    with pytest.raises(RuntimeError):
        bed.gather([("fam0", "")])


def test_snp_luts_round_trip(snplut: str):
    loaded = gelexy.load_snp_luts(snplut)
    assert list(loaded) == [gelexy.GeneticMode.A]
    np.testing.assert_allclose(loaded[gelexy.GeneticMode.A], centred_luts(DOSAGES))


def test_design_from_luts_computes_gebv(bfile: str, snplut: str):
    bed = gelexy.open_bed(bfile)
    design = gelexy.make_genetic_design(bed, gelexy.load_snp_luts(snplut))
    assert design.num_samples == 5
    assert design.num_markers == 4
    assert design.modes == [gelexy.GeneticMode.A]
    assert design.contains(gelexy.GeneticMode.A)
    assert not design.contains(gelexy.GeneticMode.D)
    np.testing.assert_allclose(design.a1_frequency[2], 0.5)

    projection = design.projection(gelexy.GeneticMode.A)
    np.testing.assert_allclose(projection.snp_luts, centred_luts(DOSAGES))

    coefficients = np.array([0.5, -1.0, 2.0, 0.25])
    expected = centred_design(DOSAGES) @ coefficients
    np.testing.assert_allclose(gelexy.gebv(projection, coefficients), expected)
    np.testing.assert_allclose(
        gelexy.gebv(projection, coefficients.astype(np.float32)), expected, rtol=1e-6
    )

    with pytest.raises(RuntimeError):
        gelexy.gebv(projection, coefficients[:3])
    with pytest.raises(RuntimeError):
        design.projection(gelexy.GeneticMode.D)


def test_gathered_bed_predicts_a_subset(bfile: str, snplut: str):
    luts = gelexy.load_snp_luts(snplut)
    bed = gelexy.open_bed(bfile)
    full = gelexy.gebv(
        gelexy.make_genetic_design(bed, luts).projection(gelexy.GeneticMode.A),
        np.array([1.0, -0.5, 0.0, 2.0]),
    )

    order = [3, 0, 2]
    bed.gather([bed.sample_ids[i] for i in order])
    subset = gelexy.gebv(
        gelexy.make_genetic_design(bed, luts).projection(gelexy.GeneticMode.A),
        np.array([1.0, -0.5, 0.0, 2.0]),
    )
    np.testing.assert_allclose(subset, full[order])

    with pytest.raises(RuntimeError):
        gelexy.make_genetic_design(bed, {})
    with pytest.raises(RuntimeError):
        gelexy.make_genetic_design(bed, {gelexy.GeneticMode.A: np.zeros((4, 3))})


def test_projection_keeps_design_alive(bfile: str, snplut: str):
    projection = gelexy.make_genetic_design(
        gelexy.open_bed(bfile), gelexy.load_snp_luts(snplut)
    ).projection(gelexy.GeneticMode.A)
    gc.collect()

    values = gelexy.gebv(projection, np.ones(4))
    assert values.shape == (5,)
    assert np.all(np.isfinite(values))
