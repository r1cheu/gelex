#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <armadillo>
#include <type_traits>

namespace nb = nanobind;

NAMESPACE_BEGIN(NB_NAMESPACE)
NAMESPACE_BEGIN(detail)

template <typename T>
constexpr bool is_mat_v = arma::is_Mat_only<T>::value;

template <typename T>
constexpr bool is_cube_v = arma::is_Cube<T>::value;

template <typename T>
constexpr int ndim_v = is_cube_v<T> ? 3 : (is_mat_v<T> ? 2 : 1);

template <typename T, typename eT = typename T::elem_type>
constexpr bool is_arma_dense_v = std::is_base_of_v<arma::Mat<eT>, T>;

template <typename T, typename eT = typename T::elem_type>
constexpr bool is_arma_sparse_v = std::is_base_of_v<arma::SpMat<eT>, T>;

// force to f-contig, since the c-contig not supported by Armadillo
template <typename T, typename eT = typename T::elem_type>
using array_for_arma_t = nb::ndarray<
    std::conditional_t<std::is_const_v<T>, const eT, eT>,
    nb::numpy,
    nb::ndim<ndim_v<T>>,
    nb::f_contig,
    nb::device::cpu>;

template <typename T>
struct type_caster<
    T,
    enable_if_t<
        is_arma_dense_v<T> && is_ndarray_scalar_v<typename T::elem_type>>>
{
    using eT = typename T::elem_type;
    using NDArray = array_for_arma_t<T>;
    using NDArrayCaster = make_caster<NDArray>;

    NB_TYPE_CASTER(T, NDArrayCaster::Name)

    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept
    {
        make_caster<NDArray> caster;
        if (!caster.from_python(
                src, flags & ~(uint8_t)cast_flags::accepts_none, cleanup))
            return false;

        NDArray& t = caster.value;
        if constexpr (ndim_v<T> == 3)
        {
            value
                = T(t.data(), t.shape(0), t.shape(1), t.shape(2), false, true);
        }
        else if constexpr (ndim_v<T> == 2)
        {
            value = T(t.data(), t.shape(0), t.shape(1), false, true);
        }
        else
        {
            value = T(t.data(), t.shape(0), false, true);
        }
        return true;
    }

    template <typename T2>
    static handle
    from_cpp(T2&& v, rv_policy policy, cleanup_list* cleanup) noexcept
    {
        policy = infer_policy<T2>(policy);
        if constexpr (std::is_pointer_v<T2>)
        {
            return from_cpp_internal((const T&)*v, policy, cleanup);
        }
        else
        {
            return from_cpp_internal((const T&)v, policy, cleanup);
        }
    }

    static handle from_cpp_internal(
        const T& v,
        rv_policy policy,
        cleanup_list* cleanup) noexcept
    {
        size_t shape[ndim_v<T>];
        int64_t strides[ndim_v<T>];

        if constexpr (ndim_v<T> == 3)
        {
            shape[0] = v.n_rows;
            shape[1] = v.n_cols;
            shape[2] = v.n_slices;
            strides[0] = 1;
            strides[1] = v.n_rows;
            strides[2] = v.n_rows * v.n_cols;
        }
        else if constexpr (ndim_v<T> == 2)
        {
            shape[0] = v.n_rows;
            shape[1] = v.n_cols;
            strides[0] = 1;
            strides[1] = v.n_rows;
        }
        else
        {
            shape[0] = v.n_elem;
            strides[0] = 1;
        }

        void* ptr = (void*)v.memptr();

        if (policy == rv_policy::move)
        {
            if ((size_t)v.n_alloc < (1024 / sizeof(eT)))
                policy = rv_policy::copy;
        }

        object owner;
        if (policy == rv_policy::move)
        {
            T* temp = new T(std::move(v));
            owner = capsule(temp, [](void* p) noexcept { delete (T*)p; });
            ptr = temp->memptr();
            policy = rv_policy::reference;
        }
        else if (policy == rv_policy::reference_internal && cleanup->self())
        {
            owner = borrow(cleanup->self());
            policy = rv_policy::reference;
        }

        object o = steal(
            NDArrayCaster::from_cpp(
                NDArray(ptr, ndim_v<T>, shape, owner, strides),
                policy,
                cleanup));

        return o.release();
    }
};

NAMESPACE_END(detail)
NAMESPACE_END(NB_NAMESPACE)
