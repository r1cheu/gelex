
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <carma>
#include "chenx/dataset/encode.h"
#include "chenx/dataset/grm.h"
#include "pybind11/iostream.h"
#include "reml.h"

namespace binds
{
using namespace arma;
using namespace carma;
using namespace chenx;
namespace py = pybind11;

py::array_t<double> hybird_value(
    const py::array_t<double>& genotype_arr,
    const py::array_t<double>& phenotype_arr)
{
    Mat<double> genotype = carma::arr_to_mat_view(genotype_arr);
    Col<double> phenotype = carma::arr_to_col_view(phenotype_arr);
    return carma::mat_to_arr(chenx::hybird_value(genotype, phenotype));
}

void hybird(
    py::array_t<double>& genotype_arr,
    const py::array_t<double>& hybird_value_arr)
{
    Mat<double> genotype = carma::arr_to_mat(genotype_arr);
    Mat<double> hybird_value = carma::arr_to_mat_view(hybird_value_arr);
    chenx::hybird(genotype, hybird_value);
}

py::array_t<double> Amat(py::array_t<double>& genotype_arr)
{
    Mat<double> genotype_mat = carma::arr_to_mat(genotype_arr);
    return carma::mat_to_arr(chenx::Amat(genotype_mat));
}

py::array_t<double> Dmat(py::array_t<double>& genotype_arr)
{
    Mat<double> genotype_mat = carma::arr_to_mat(genotype_arr);
    return carma::mat_to_arr(chenx::Dmat(genotype_mat));
}

py::array_t<double> Amat_rbf(const py::array_t<double>& X_arr, double bandwidth)
{
    Mat<double> X = carma::arr_to_mat_view(X_arr);
    return carma::mat_to_arr(Add_rbf_kernel(X, bandwidth));
}

py::array_t<double> Dmat_rbf(py::array_t<double>& X_arr, double bandwidth)
{
    Mat<double> X = carma::arr_to_mat(X_arr);
    return carma::mat_to_arr(Dom_rbf_kernel(X, bandwidth));
}

void bind_dataset(py::module& m)
{
    m.def(
        "_hybrid_value",
        &hybird_value,
        "calculate hybird value",
        py::arg("genotype"),
        py::arg("phenotype"));
    m.def(
        "_hybrid",
        &hybird,
        "Replace D value",
        py::arg("genotype"),
        py::arg("hybird_value"));
    m.def(
        "_Amat",
        &Amat,
        "Additive Genetic Relationship Matrix",
        py::arg("genotype"));
    m.def(
        "_Dmat",
        &Dmat,
        "Dominance Genetic Relationship Matrix",
        py::arg("genotype"));
    m.def(
        "_Amat_rbf",
        &Amat_rbf,
        "Additive Genetic Relationship Matrix with RBF kernel",
        py::arg("X"),
        py::arg("bandwidth"));
    m.def(
        "_Dmat_rbf",
        &Dmat_rbf,
        "Dominance Genetic Relationship Matrix with RBF kernel",
        py::arg("X"),
        py::arg("bandwidth"));
}

PYBIND11_MODULE(_core, m)
{
    bind_dataset(m);
    // bind Additive Genetic Relationship Matrix
    py::class_<REMLLoop<double>>(m, "REMLLoop")
        .def(py::init<
             const py::array_t<double>&,
             const py::array_t<double>&,
             const py::array_t<uword>&,
             const py::array_t<double>&,
             std::vector<std::string>>())
        .def(
            "run",
            [](REMLLoop<double>& self,
               const py::array_t<double>& varcomp_prior_arr,
               std::string_view method,
               bool em_init,
               size_t max_iteration,
               double tolerance,
               bool verbose)
            {
                // Redirect std::cout to sys.stdout
                py::scoped_ostream_redirect redirect(
                    std::cout, py::module_::import("sys").attr("stdout"));
                self.run(
                    varcomp_prior_arr,
                    method,
                    em_init,
                    max_iteration,
                    tolerance,
                    verbose);
            },
            py::arg("varcomp"),
            py::arg("method"),
            py::arg("em_init"),
            py::arg("max_iteration"),
            py::arg("tolerance"),
            py::arg("verbose"))
        .def("get_varcomp", &REMLLoop<double>::get_varcomp)
        .def("get_beta", &REMLLoop<double>::get_beta)
        .def("get_blup", &REMLLoop<double>::get_blup)
        .def("get_gebv", &REMLLoop<double>::get_gebv, py::arg("full_X_arr"));
}
}  // namespace binds
