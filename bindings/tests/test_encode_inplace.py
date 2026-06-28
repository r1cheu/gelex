import gelex
import numpy as np
import pytest


def make_genotypes():
    # rows = samples, cols = markers; raw genotypes 0/1/2
    return np.asfortranarray(
        np.array([[0, 1, 2], [2, 1, 0], [1, 0, 2]], dtype=np.float64)
    )


def test_enums_exposed():
    assert gelex.GeneticMode.A != gelex.GeneticMode.D
    assert gelex.GenotypeMethod.Standardize is not None


def test_encode_inplace_writes_correct_values_in_place():
    # End-to-end oracle: additive Center subtracts each column's dosage mean.
    # Computing the expectation independently in NumPy proves the binding wires
    # the args, the double instantiation, and the in-place writeback correctly
    # -- not merely that *some* mutation happened. Exhaustive per-method numeric
    # coverage lives in tests/test_locus_encoding.cpp.
    g = np.asfortranarray(np.array([[0.0], [1.0], [2.0]], dtype=np.float64))

    gelex.encode_inplace(g, gelex.GeneticMode.A, gelex.GenotypeMethod.Center)

    assert np.allclose(g.ravel(), [-1.0, 0.0, 1.0])  # column mean (1.0) removed


def test_return_object_structure_readable():
    # the nested return type is exposed: spec, loci list, std::array code,
    # and the nested LocusStats are all reachable from Python.
    g = make_genotypes()
    enc = gelex.encode_inplace(g, gelex.GeneticMode.A, gelex.GenotypeMethod.Standardize)

    assert isinstance(enc, gelex.LociEncoding)
    assert enc.spec.effect == gelex.GeneticMode.A
    assert len(enc.loci) == g.shape[1]
    assert len(enc.loci[0].code) == 3
    assert isinstance(enc.loci[0].stats, gelex.LocusStats)


def test_marker_offset_applied():
    g = make_genotypes()
    enc = gelex.encode_inplace(
        g,
        gelex.GeneticMode.A,
        gelex.GenotypeMethod.Center,
        marker_offset=10,
    )
    assert enc.loci[0].marker_index == 10


def test_c_order_rejected():
    g = np.ascontiguousarray(np.array([[0, 1, 2], [2, 1, 0]], dtype=np.float64))
    with pytest.raises(TypeError):
        gelex.encode_inplace(g, gelex.GeneticMode.A, gelex.GenotypeMethod.Standardize)


def test_wrong_dtype_rejected():
    g = np.asfortranarray(np.array([[0, 1, 2], [2, 1, 0]], dtype=np.float32))
    with pytest.raises(TypeError):
        gelex.encode_inplace(g, gelex.GeneticMode.A, gelex.GenotypeMethod.Standardize)
