import gelex
import numpy as np
import pytest


def test_enums_exposed():
    assert gelex.GeneticMode.A != gelex.GeneticMode.D
    assert gelex.GenotypeMethod.Standardize is not None


def test_encode_inplace_uses_gelex_encoding():
    genotypes = np.asfortranarray(
        np.array([[0.0], [1.0], [2.0]], dtype=np.float64)
    )

    result = gelex.encode_inplace(
        genotypes,
        gelex.GeneticMode.A,
        gelex.GenotypeMethod.Center,
    )

    assert result is None
    np.testing.assert_allclose(genotypes[:, 0], [-1.0, 0.0, 1.0])


def test_encode_inplace_rejects_incompatible_arrays():
    c_order = np.array([[0.0, 1.0], [1.0, 2.0]], order="C")
    float32 = np.array([[0.0, 1.0], [1.0, 2.0]], dtype=np.float32, order="F")
    readonly = np.array([[0.0, 1.0], [1.0, 2.0]], order="F")
    readonly.flags.writeable = False

    for genotypes in (c_order, float32, readonly):
        with pytest.raises(TypeError):
            gelex.encode_inplace(
                genotypes,
                gelex.GeneticMode.A,
                gelex.GenotypeMethod.Center,
            )
