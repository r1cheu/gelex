from pathlib import Path

import arviz as az
import gelex
import numpy as np
import pytest

DRAWS = 20

# Mirrors an MCMC run with 2 fixed columns, one random block of 3 levels,
# 5 markers and 20 draws; cell (row, draw) = row + draw / 100.
LAYOUT = [
    ("fixed/coefficients", gelex.BinaryType.float32, 2),
    ("random/batch/coefficients", gelex.BinaryType.float32, 3),
    ("random/batch/variance", gelex.BinaryType.float64, 1),
    ("genetic/A/coefficients", gelex.BinaryType.float32, 5),
    ("genetic/A/variance", gelex.BinaryType.float64, 1),
    ("genetic/A/assignment", gelex.BinaryType.uint8, 5),
    ("genetic/A/probability", gelex.BinaryType.float64, 1),
    ("residual/variance", gelex.BinaryType.float64, 1),
    ("genetic/A/explained_variance", gelex.BinaryType.float64, 1),
    ("genetic/A/heritability", gelex.BinaryType.float64, 1),
    ("genetic/total/explained_variance", gelex.BinaryType.float64, 1),
    ("genetic/total/heritability", gelex.BinaryType.float64, 1),
]

NUMPY_DTYPE = {
    gelex.BinaryType.float64: np.float64,
    gelex.BinaryType.float32: np.float32,
    gelex.BinaryType.uint8: np.uint8,
}


def expected(rows: int, dtype) -> np.ndarray:
    """(rows, DRAWS) in column-major order, matching the on-disk layout."""
    grid = np.arange(rows)[:, None] + np.arange(DRAWS)[None, :] / 100.0
    return np.asfortranarray(grid.astype(dtype))


@pytest.fixture
def draws_path(tmp_path: Path) -> Path:
    path = tmp_path / "fixture.draws"
    with gelex.DenseWriter(str(path)) as writer:
        payloads = [
            (
                writer.reserve(name, dtype, (rows, DRAWS)),
                expected(rows, NUMPY_DTYPE[dtype]),
            )
            for name, dtype, rows in LAYOUT
        ]
        for draw in range(DRAWS):
            for payload, values in payloads:
                payload.append(values[:, draw])
    return path


def test_writer_round_trips_through_reader(tmp_path: Path):
    path = tmp_path / "roundtrip.draws"
    writer = gelex.DenseWriter(str(path))
    appended = writer.reserve("appended", gelex.BinaryType.float32, (2, 3))
    whole = writer.reserve("whole", gelex.BinaryType.float64, (2, 3))
    assert appended.identifier == "appended"

    for column in np.array([[1, 2], [3, 4], [5, 6]], dtype=np.float32):
        appended.append(column)
    matrix = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    with pytest.raises(TypeError):
        whole.write(matrix)
    whole.write(np.asfortranarray(matrix))
    assert writer.is_open
    writer.close()
    assert not writer.is_open

    reader = gelex.BinaryReader(str(path))
    np.testing.assert_array_equal(
        reader["appended"], np.array([[1, 3, 5], [2, 4, 6]], dtype=np.float32)
    )
    np.testing.assert_array_equal(reader["whole"], matrix)


def test_writer_validates_columns(tmp_path: Path):
    writer = gelex.DenseWriter(str(tmp_path / "invalid.draws"))
    payload = writer.reserve("p", gelex.BinaryType.float64, (2, 1))

    assert isinstance(payload, gelex.DenseStreamF64)
    with pytest.raises(TypeError):
        payload.append(np.zeros(2, dtype=np.float32))
    with pytest.raises(Exception, match="expected 2 column values"):
        payload.append(np.zeros(3))
    with pytest.raises(Exception, match="duplicate|already"):
        writer.reserve("p", gelex.BinaryType.float64, (1, 1))

    payload.append(np.zeros(2))
    writer.close()
    with pytest.raises(Exception, match="closed"):
        payload.append(np.zeros(2))
    with pytest.raises(Exception, match="closed"):
        writer.reserve("q", gelex.BinaryType.float64, (1, 1))


def test_reader_lists_payloads(draws_path: Path):
    reader = gelex.BinaryReader(str(draws_path))

    assert len(reader) == len(LAYOUT)
    assert set(reader.keys()) == {name for name, _, _ in LAYOUT}
    assert [info.identifier for info in reader.payloads()] == reader.keys()
    assert "residual/variance" in reader
    assert "missing" not in reader

    info = reader.info("genetic/A/assignment")
    assert info.identifier == "genetic/A/assignment"
    assert info.type == gelex.BinaryType.uint8
    assert info.shape == (5, DRAWS)


def test_reader_exposes_readonly_column_major_views(draws_path: Path):
    reader = gelex.BinaryReader(str(draws_path))

    for name, dtype, rows in LAYOUT:
        view = reader[name]
        assert view.dtype == NUMPY_DTYPE[dtype]
        assert view.shape == (rows, DRAWS)
        assert view.flags.f_contiguous
        assert not view.flags.writeable
        np.testing.assert_array_equal(view, expected(rows, NUMPY_DTYPE[dtype]))

    with pytest.raises(Exception):
        reader["missing"]


def test_views_keep_the_reader_alive(draws_path: Path):
    view = gelex.BinaryReader(str(draws_path))["fixed/coefficients"]
    np.testing.assert_array_equal(view, expected(2, np.float32))


def test_read_draws_skips_marker_sized_payloads_by_default(draws_path: Path):
    idata = gelex.read_draws(draws_path)
    posterior = idata.posterior

    assert "genetic.A.coefficients" not in posterior
    assert "genetic.A.assignment" not in posterior
    assert posterior["residual.variance"].shape == (1, DRAWS)
    assert posterior["random.batch.coefficients"].shape == (1, DRAWS, 3)
    np.testing.assert_allclose(
        posterior["random.batch.coefficients"].values[0],
        expected(3, np.float32).T,
    )

    summary = az.summary(idata, var_names=["residual.variance"])
    assert "ess_bulk" in summary.columns


def test_read_draws_includes_markers_on_request(draws_path: Path):
    posterior = gelex.read_draws(draws_path, include_markers=True).posterior

    assert posterior["genetic.A.coefficients"].shape == (1, DRAWS, 5)
    assert posterior["genetic.A.assignment"].shape == (1, DRAWS, 5)
