#pragma once

#include <nanobind/ndarray.h>
#include "dense.h"

#include <memory>
#include <type_traits>
#include <utility>

NAMESPACE_BEGIN(NB_NAMESPACE)
NAMESPACE_BEGIN(detail)

template <typename T>
struct type_caster<T, enable_if_t<is_arma_sparse_v<T>>>
{
    using eT = typename T::elem_type;
    using Index = arma::uword;

    using NDArray = ndarray<numpy, eT, shape<-1>>;
    using IndexNDArray = ndarray<numpy, Index, shape<-1>>;

    using NDArrayCaster = make_caster<NDArray>;
    using IndexCaster = make_caster<IndexNDArray>;

    NB_TYPE_CASTER(
        T,
        const_name("scipy.sparse.csc_matrix[") + make_caster<eT>::Name
            + const_name("]"))

    NDArrayCaster data_caster;
    IndexCaster indices_caster, indptr_caster;

    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept
    {
        object obj = borrow(src);

        try
        {
            object matrix_type
                = module_::import_("scipy.sparse").attr("csc_matrix");
            if (!obj.type().is(matrix_type))
                obj = matrix_type(obj);

            if (!cast<bool>(obj.attr("has_sorted_indices")))
                obj.attr("sort_indices")();

            if (object data_o = obj.attr("data");
                !data_caster.from_python(data_o, flags, cleanup))
                return false;

            if (object indices_o = obj.attr("indices");
                !indices_caster.from_python(indices_o, flags, cleanup))
                return false;

            if (object indptr_o = obj.attr("indptr");
                !indptr_caster.from_python(indptr_o, flags, cleanup))
                return false;

            object shape_o = obj.attr("shape");
            if (len(shape_o) != 2)
                return false;

            Index rows = cast<Index>(shape_o[0]),
                  cols = cast<Index>(shape_o[1]),
                  nnz = cast<Index>(obj.attr("nnz"));

            value = arma::SpMat<eT>(
                arma::uvec(indices_caster.value.data(), nnz),
                arma::uvec(indptr_caster.value.data(), cols + 1),
                arma::Col<eT>(data_caster.value.data(), nnz),
                rows,
                cols,
                false);
            return true;
        }
        catch (const python_error&)
        {
            return false;
        }
    }

    static handle
    from_cpp(T&& v, rv_policy policy, cleanup_list* cleanup) noexcept
    {
        if (policy == rv_policy::automatic
            || policy == rv_policy::automatic_reference)
            policy = rv_policy::move;

        return from_cpp((const T&)v, policy, cleanup);
    }

    template <typename T2>
    static handle
    from_cpp(T2&& v, rv_policy policy, cleanup_list* cleanup) noexcept
    {
        policy = infer_policy<T2>(policy);
        if constexpr (std::is_pointer_v<T2>)
            return from_cpp_internal((const T&)*v, policy, cleanup);
        else
            return from_cpp_internal((const T&)v, policy, cleanup);
    }

    static handle
    from_cpp_internal(const T& v, rv_policy policy, cleanup_list*) noexcept
    {
        object matrix_type;
        try
        {
            matrix_type = module_::import_("scipy.sparse").attr("csc_matrix");
        }
        catch (python_error& e)
        {
            e.restore();
            return handle();
        }

        const Index rows = v.n_rows, cols = v.n_cols;
        const size_t data_shape[] = {(size_t)v.n_nonzero};
        const size_t indices_shape[] = {(size_t)(cols + 1)};

        T* src = std::addressof(const_cast<T&>(v));
        object owner;
        if (policy == rv_policy::move)
        {
            src = new T(std::move(v));
            owner = capsule(src, [](void* p) noexcept { delete (T*)p; });
        }

        NDArray data(const_cast<eT*>(src->values), 1, data_shape, owner);
        IndexNDArray indices(
            const_cast<Index*>(src->row_indices), 1, data_shape, owner);
        IndexNDArray indptr(
            const_cast<Index*>(src->col_ptrs), 1, indices_shape, owner);
        try
        {
            return matrix_type(
                       nanobind::make_tuple(
                           std::move(data),
                           std::move(indices),
                           std::move(indptr)),
                       nanobind::make_tuple(rows, cols))
                .release();
        }
        catch (python_error& e)
        {
            e.restore();
            return handle();
        }
    }
};

NAMESPACE_END(detail)
NAMESPACE_END(NB_NAMESPACE)
