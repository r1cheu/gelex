#pragma once
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <armadillo>
#include <cassert>
#include <type_traits>
#include <utility>

namespace bind
{
namespace nb = nanobind;
template <typename T, int NDIM, typename CM = nb::any_contig>
using ndarray_type
    = nb::ndarray<T, nb::numpy, nb::ndim<NDIM>, nb::device::cpu, CM>;

template <typename T>
using vec = ndarray_type<T, 1>;
template <typename T>
using vec_view = ndarray_type<const T, 1>;
template <typename T>
using mat = ndarray_type<T, 2, nb::f_contig>;
template <typename T>
using mat_view = ndarray_type<const T, 2, nb::f_contig>;
template <typename T>
using cube = ndarray_type<T, 3, nb::f_contig>;
template <typename T>
using cube_view = ndarray_type<const T, 3, nb::f_contig>;

using dmat = mat<double>;
using dvec = vec<double>;
using dcube = cube<double>;

template <typename ArmaType>
struct PyOwned;
template <typename ArmaType>
struct PyView;
template <typename PyType>
struct Arma;

template <typename T>
struct PyOwned<arma::Col<T>>
{
    using type = vec<T>;
    static constexpr int ndim = 1;
};

template <typename T>
struct PyView<arma::Col<T>>
{
    using type = vec_view<T>;
    static constexpr int ndim = 1;
};
template <typename T>
struct PyOwned<arma::Mat<T>>
{
    using type = mat<T>;
    static constexpr int ndim = 2;
};
template <typename T>
struct PyView<arma::Mat<T>>
{
    using type = mat_view<T>;
    static constexpr int ndim = 2;
};
template <typename T>
struct PyOwned<arma::Cube<T>>
{
    using type = cube<T>;
    static constexpr int ndim = 3;
};
template <typename T>
struct PyView<arma::Cube<T>>
{
    using type = cube_view<T>;
    static constexpr int ndim = 3;
};

template <typename T>
struct Arma<vec<T>>
{
    using type = arma::Col<T>;
    static constexpr int ndim = 1;
};
template <typename T>
struct Arma<mat<T>>
{
    using type = arma::Mat<T>;
    static constexpr int ndim = 2;
};
template <typename T>
struct Arma<cube<T>>
{
    using type = arma::Cube<T>;
    static constexpr int ndim = 3;
};

template <typename ArmaType>
static void Deleter(void* p) noexcept
{
    delete static_cast<ArmaType*>(p);
}

template <typename Array, int NDIM, typename ArmaObj>
Array ToPyImpl(ArmaObj* arma_ptr, nb::capsule* deleter = nullptr)
{
    if constexpr (NDIM == 1)
    {
        return deleter ? Array(arma_ptr->memptr(), {arma_ptr->n_elem}, *deleter)
                       : Array(arma_ptr->memptr(), {arma_ptr->n_elem});
    }
    else if constexpr (NDIM == 2)
    {
        return deleter ? Array(
                             arma_ptr->memptr(),
                             {arma_ptr->n_rows, arma_ptr->n_cols},
                             *deleter)
                       : Array(
                             arma_ptr->memptr(),
                             {arma_ptr->n_rows, arma_ptr->n_cols});
    }
    else
    {
        return deleter ? Array(
                             arma_ptr->memptr(),
                             {arma_ptr->n_rows,
                              arma_ptr->n_cols,
                              arma_ptr->n_slices},
                             *deleter)
                       : Array(
                             arma_ptr->memptr(),
                             {arma_ptr->n_rows,
                              arma_ptr->n_cols,
                              arma_ptr->n_slices});
    }
}

template <typename ArmaType>
auto ToPy(ArmaType&& arma_obj)
{
    using DecayType = std::decay_t<ArmaType>;
    if constexpr (std::is_rvalue_reference_v<ArmaType&&>)
    {
        using ArrayType = typename PyOwned<DecayType>::type;
        constexpr int ndim = PyOwned<DecayType>::ndim;
        auto* ptr = new DecayType(std::forward<ArmaType>(arma_obj));
        nb::capsule capsule(ptr, Deleter<DecayType>);
        return ToPyImpl<ArrayType, ndim>(ptr, &capsule);
    }
    else if constexpr (
        std::is_lvalue_reference_v<ArmaType&&> || std::is_const_v<ArmaType&&>)
    {
        using ArrayType = typename PyView<DecayType>::type;
        constexpr int ndim = PyView<DecayType>::ndim;
        return ToPyImpl<ArrayType, ndim>(&arma_obj);
    }
}

template <typename PyType>
auto ToArma(PyType array)
{
    using DecayType = std::decay_t<PyType>;
    using ArmaType = typename Arma<DecayType>::type;
    constexpr int ndim = Arma<DecayType>::ndim;

    if constexpr (ndim == 1)
    {
        return ArmaType(array.data(), array.shape(0), false, true);
    }
    else if constexpr (ndim == 2)
    {
        return ArmaType(
            array.data(), array.shape(0), array.shape(1), false, true);
    }
    else
    {
        return ArmaType(
            array.data(),
            array.shape(0),
            array.shape(1),
            array.shape(2),
            false,
            true);
    }
}
}  // namespace bind
