import gc

import gelexy._test as t
import numpy as np
import pytest
from numpy.testing import assert_array_equal
from scipy import sparse


def test01_vec():
    a = np.array([1, 2, 3], dtype=np.int64)
    b = np.array([0, 1, 2], dtype=np.int64)
    c = np.array([1, 3, 5], dtype=np.int64)
    x = np.arange(10000, dtype=np.int64)
    af = np.float32(a)
    bf = np.float32(b)

    # Check call with dynamically sized arrays
    assert_array_equal(t.add_ivec(a, b), c)

    # Implicit conversion
    assert_array_equal(t.add_ivec(af, b), c)

    with pytest.raises(TypeError, match="incompatible function arguments"):
        t.add_ivec(a, bf)

    # Try with a big array. This will move the result to avoid a copy
    assert_array_equal(t.add_ivec(x, x), 2 * x)


def test02_mat():
    a = np.array([[1, 2], [3, 4]], dtype=np.int64, order="F")
    b = np.array([[0, 1], [2, 3]], dtype=np.int64, order="F")
    c = np.array([[1, 3], [5, 7]], dtype=np.int64, order="F")
    x = np.arange(10000).reshape(100, 100).astype(np.int64, order="F")

    af = np.float32(a)
    bf = np.float32(b)
    ac = np.array([[1, 2], [3, 4]], dtype=np.int64, order="C")
    bc = np.array([[0, 1], [2, 3]], dtype=np.int64, order="C")

    assert_array_equal(t.add_imat(a, b), c)

    # Implicit conversion type
    assert_array_equal(t.add_imat(af, b), c)
    # implicit conversion order
    assert_array_equal(t.add_imat(ac, b), c)

    with pytest.raises(TypeError, match="incompatible function arguments"):
        t.add_imat(a, bf)
        t.add_imat(a, bc)

    # Try with a big array. This will move the result to avoid a copy
    assert_array_equal(t.add_imat(x, x), 2 * x)


def test03_inplace_update():
    a = np.array([1, 2, 3], dtype=np.int64)
    b = np.array([[0, 1], [2, 3]], dtype=np.int64, order="F")

    c = np.array([99, 2, 3], dtype=np.int64)
    d = np.array([[99, 1], [2, 3]], dtype=np.int64, order="F")

    t.update_ivec(a)
    assert_array_equal(a, c)

    t.update_imat(b)
    assert_array_equal(b, d)


def test04_prop():
    for j in range(3):
        c = t.ClassWithArmaMember()
        ref = np.ones((2, 2))
        if j == 0:
            c.member = ref

        for i in range(2):
            member = c.member
            if j == 2 and i == 0:
                member[0, 0] = 10
                ref[0, 0] = 10
            assert_array_equal(member, ref)
            del member
            gc.collect()
            gc.collect()

        member = c.member
        assert_array_equal(c.member_ro_ref, ref)
        assert_array_equal(c.member_ro_copy, ref)
        del c
        gc.collect()
        gc.collect()
        assert_array_equal(member, ref)


def test05_init_assert():
    a = np.array([[1, 2, 3], [1, 2, 3]], dtype=np.float64, order="F")
    b = np.array([[0, 2, 3], [1, 2, 3]], dtype=np.float64, order="F")

    A = t.ClassInitFromPython(a)
    A.get_mat()[0, 0] = 0

    # assert a to make sure use same address
    assert_array_equal(a, b)

    del a
    gc.collect()
    gc.collect()
    assert_array_equal(A.get_mat(), b)


def test06_optional_class():
    a = np.array([[1, 2, 3], [1, 2, 3]], dtype=np.float64, order="F")
    b = np.array([[0, 2, 3], [1, 2, 3]], dtype=np.float64, order="F")

    A = t.ClassInitOptional(a)
    A.get_mat()[0, 0] = 0
    # assert a to make sure use same address
    assert_array_equal(a, b)

    A = t.ClassInitOptional()
    assert_array_equal(A.get_mat(), b)


def test07_sparse():
    a = np.array([[1, 3, 4], [0, 0, 0], [1, 0, 0]], dtype=np.float64, order="F")
    assert_array_equal(a, t.sparse().todense())

    b = a + 1
    a = sparse.csc_matrix(a)
    assert_array_equal(b, t.sparse_add(a))
