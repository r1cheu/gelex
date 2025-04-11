import gelexy._test as t
import numpy as np
import pytest


def test_vector_ptr_consistency():
    np_array = np.array([1.0, 2.0, 3.0], dtype=np.float64)

    py_ptr = np_array.__array_interface__["data"][0]
    cpp_ptr = t.get_vec_ptr_value(np_array)

    assert py_ptr == cpp_ptr, (
        f"Vector pointer mismatch: Python={py_ptr}, C++={cpp_ptr}"
    )

    np_array[0] = 42.0
    new_array = np.array([42.0, 2.0, 3.0], dtype=np.float64)
    assert np.array_equal(np_array, new_array)


def test_matrix_ptr_consistency():
    np_mat = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64, order="F")
    py_ptr = np_mat.__array_interface__["data"][0]
    cpp_ptr = t.get_mat_ptr_value(np_mat)

    assert py_ptr == cpp_ptr, (
        f"Matrix pointer mismatch: Python={py_ptr}, C++={cpp_ptr}"
    )

    np_mat[0, 0] = 42.0
    modified = np.array([[42.0, 2.0], [3.0, 4.0]], dtype=np.float64, order="F")
    assert np.array_equal(np_mat, modified)


def test_optional_matrix():
    np_mat = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64, order="F")
    py_ptr = np_mat.__array_interface__["data"][0]

    cpp_ptr = t.get_optional_mat_ptr_value(np_mat)
    assert py_ptr == cpp_ptr, (
        f"Optional matrix pointer mismatch: Python={py_ptr}, C++={cpp_ptr}"
    )

    none_ptr = t.get_optional_mat_ptr_value(None)
    assert none_ptr is None, f"Expected None, got {none_ptr}"


def test_return_mat_lvalue():
    a = t.return_A()
    result_mat = a.const_arr()
    assert result_mat[0, 0] == 99.0

    result_addr = result_mat.__array_interface__["data"][0]
    assert result_addr == a.get_arr_ptr(), (
        f"Lvalue reference address mismatch: cpp={a.get_arr_ptr()}, Returned={result_addr}"
    )

    with pytest.raises(ValueError, match="assignment destination is read-only"):
        result_mat[0, 0] = 55.0

    result_mat2 = a.mutable_arr()
    result_mat2[0, 0] = 77
    assert a.const_arr()[0, 0] == 77.0, (
        "Change of reference wasn't reflected in the original matrix"
    )


def test_return_mat_rvalue():
    result_mat = t.return_mat_rvalue()
    assert result_mat[0, 0] == 99.0
    result_mat[0, 0] = 55.0
    result_mat2 = t.return_mat_rvalue()
    assert result_mat2[0, 0] == 99.0, "Rvalue wasn't properly moved/copied"


def test_return_rvalue():
    a = t.return_A()
    result = a.new_arr()
    result_addr = result.__array_interface__["data"][0]
    cpp_addr = a.new_arr_ptr()
    assert result_addr == cpp_addr, (
        "Rvalue address mismatch: cpp={cpp_addr}, Returned={result_addr}"
    )


def test_matrix_data_integrity():
    value = 3.14
    test_mat = t.create_test_matrix(value)

    assert np.allclose(test_mat, np.ones((3, 3)) * value)

    test_mat[0, 0] = 2.71
    assert test_mat[0, 0] == 2.71


def test_modify_and_return():
    initial_mat = np.ones((3, 3), dtype=np.float64, order="F")

    initial_addr = initial_mat.__array_interface__["data"][0]
    result_mat = t.modify_and_return_matrix(initial_mat)

    result_addr = result_mat.__array_interface__["data"][0]
    assert initial_addr == result_addr, (
        "Returned matrix should be the same object"
    )
    assert np.allclose(result_mat, np.ones((3, 3)) * 2.0)
    assert np.allclose(initial_mat, np.ones((3, 3)) * 2.0)
