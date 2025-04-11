#pragma once
#include <cassert>
#include <type_traits>
#include <utility>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <armadillo>

namespace bind
{
namespace nb = nanobind;

template <typename ArmaType>
constexpr int ndim_v
    = arma::is_Row<ArmaType>::value || arma::is_Col<ArmaType>::value ? 1
      : arma::is_Mat<ArmaType>::value                                ? 2
      : arma::is_Cube<ArmaType>::value                               ? 3
                                                                     : -1;

template <typename ArmaType, typename Scalar = typename ArmaType::elem_type>
using array_for_arma_t = nb::ndarray<
    std::conditional_t<std::is_const_v<ArmaType>, const Scalar, Scalar>,
    nb::numpy,
    nb::ndim<ndim_v<ArmaType>>,
    nb::f_contig,
    nb::device::cpu>;

using arr1d = array_for_arma_t<arma::dvec>;
using arr2d = array_for_arma_t<arma::dmat>;
using arr3d = array_for_arma_t<arma::dcube>;

template <typename Numpy, typename Scalar = typename Numpy::Scalar>
struct arma_for_array;

template <>
struct arma_for_array<arr1d>
{
    using type = arma::dvec;
};

template <>
struct arma_for_array<arr2d>
{
    using type = arma::dmat;
};

template <>
struct arma_for_array<arr3d>
{
    using type = arma::dcube;
};

template <typename Numpy>
using arma_for_array_t = typename arma_for_array<Numpy>::type;

template <typename ArmaType>
auto to_py(ArmaType&& arma_obj)
{
    using Arma = std::decay_t<ArmaType>;
    using Numpy = std::conditional_t<
        std::is_const_v<std::remove_reference_t<ArmaType&&>>,
        array_for_arma_t<const Arma>,
        array_for_arma_t<Arma>>;

    if constexpr (std::is_rvalue_reference_v<ArmaType&&>)
    {
        auto* ptr = new Arma(std::forward<ArmaType>(arma_obj));
        nb::capsule owner(
            ptr, [](void* p) noexcept { delete static_cast<Arma*>(p); });
        if constexpr (ndim_v<Arma> == 1)
        {
            return Numpy(ptr->memptr(), {ptr->n_elem}, owner);
        }
        if constexpr (ndim_v<Arma> == 2)
        {
            return Numpy(ptr->memptr(), {ptr->n_rows, ptr->n_cols}, owner);
        }
        if constexpr (ndim_v<Arma> == 3)
        {
            return Numpy(
                ptr->memptr(),
                {ptr->n_rows, ptr->n_cols, ptr->n_slices},
                owner);
        }
    }
    else if constexpr (std::is_lvalue_reference_v<ArmaType&&>)
    {
        if constexpr (ndim_v<Arma> == 1)
        {
            return Numpy(arma_obj.memptr(), {arma_obj.n_elem});
        }
        if constexpr (ndim_v<Arma> == 2)
        {
            return Numpy(arma_obj.memptr(), {arma_obj.n_rows, arma_obj.n_cols});
        }
        if constexpr (ndim_v<Arma> == 3)
        {
            return Numpy(
                arma_obj.memptr(),
                {arma_obj.n_rows, arma_obj.n_cols, arma_obj.n_slices});
        }
    }
}

template <typename Numpy, typename Scalar = typename Numpy::Scalar>
auto to_arma(Numpy arr)
{
    using ArmaType = arma_for_array_t<std::decay_t<Numpy>>;
    if constexpr (arma::is_Col<ArmaType>::value)
    {
        return ArmaType(arr.data(), arr.shape(0), false, true);
    }
    else if constexpr (arma::is_Mat<ArmaType>::value)
    {
        return ArmaType(arr.data(), arr.shape(0), arr.shape(1), false, true);
    }
    else if constexpr (arma::is_Cube<ArmaType>::value)
    {
        return ArmaType(
            arr.data(), arr.shape(0), arr.shape(1), arr.shape(2), false, true);
    }
}

template <typename Numpy, typename Scalar = typename Numpy::Scalar>
arma::Row<Scalar> ToRowVec(Numpy arr)
{
    return arma::Row<Scalar>(arr.data(), arr.shape(0), false, true);
}

}  // namespace bind
