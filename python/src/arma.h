#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <type_traits>
#include <utility>

namespace bind
{
namespace nb = nanobind;

using dvec = nb::ndarray<
    double,
    nb::numpy,
    nb::ndim<1>,
    nb::device::cpu,
    nb::any_contig>;  // 1d array is no difference in memory layout

using dvec_view = nb::ndarray<
    const double,
    nb::numpy,
    nb::ndim<1>,
    nb::device::cpu,
    nb::any_contig>;

using dmat_a = nb::
    ndarray<double, nb::numpy, nb::ndim<2>, nb::device::cpu, nb::any_contig>;

using dmat_a_view = nb::ndarray<
    const double,
    nb::numpy,
    nb::ndim<2>,
    nb::device::cpu,
    nb::any_contig>;

using dmat = nb::
    ndarray<double, nb::numpy, nb::ndim<2>, nb::device::cpu, nb::f_contig>;

using dmat_view = nb::ndarray<
    const double,
    nb::numpy,
    nb::ndim<2>,
    nb::device::cpu,
    nb::f_contig>;

using dcube = nb::
    ndarray<double, nb::numpy, nb::ndim<3>, nb::device::cpu, nb::f_contig>;

using dcube_view = nb::ndarray<
    const double,
    nb::numpy,
    nb::ndim<3>,
    nb::device::cpu,
    nb::f_contig>;

template <typename ArmaType>
static void Deleter(void* p) noexcept
{
    auto* arma_ptr = static_cast<ArmaType*>(p);
    delete arma_ptr;
}

/**
 * @brief converts an temporary armadillo matrix to a python-side

 * This function takes an Armadillo matrix (or any compatible type) and converts
 * it to a Python array using nanobind. The function handles memory management
 * by creating a capsule that will delete the Armadillo matrix when the python
 * array is no longer in use.
 *
 * @tparam PyType The type of the array to python side(e.g., dmat, dvec, etc.
 * see above)
 * @tparam ArmaType The type of the Armadillo matrix (must be an rvalue
 * reference)
 * @param arma_mat The Armadillo matrix to convert
 * @return A array containing the data from the Armadillo matrix
 *
 * @note The input `arma_mat` must be an rvalue reference (e.g., a temporary
 * object or the result of std::move). This ensures that the matrix is not
 * copied unnecessarily.
 */
template <typename PyType, typename ArmaType>
static PyType ToPy(ArmaType&& arma_mat)
{
    static_assert(
        !std::is_lvalue_reference_v<ArmaType>,
        "ToNumpy: ArmaType must be rvalue reference");
    auto* mat_ptr
        = new std::decay_t<ArmaType>(std::forward<ArmaType>(arma_mat));
    nb::capsule deleter(mat_ptr, Deleter<std::decay_t<ArmaType>>);
    return PyType(
        mat_ptr->memptr(), {mat_ptr->n_rows, mat_ptr->n_cols}, deleter);
}

template <typename PyType, typename ArmaType>
static PyType ToPyView(const ArmaType& arma_mat)
{
    return PyType(arma_mat.memptr(), {arma_mat.n_rows, arma_mat.n_cols});
}
}  // namespace bind
