#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <carma>
#include "chenx/dataset/encode.h"
#include "chenx/dataset/grm.h"
#include "chenx/dataset/impute.h"

namespace binds {
using namespace arma;
namespace py = pybind11;
py::array_t<double> grm(const py::array_t<double>& X_arr, std::string encode_method) {
    dmat X = carma::arr_to_mat_view(X_arr);
    chenx::dataset::normalize(X, encode_method);
    return carma::mat_to_arr(chenx::dataset::cal_grm(X));
}

void impute(py::array_t<double>& X_arr, std::string_view method) {
    Mat<double> X = carma::arr_to_mat(X_arr, false);
    if (method == "mean") {
        chenx::dataset::mean_impute(X);
    } else if (method == "median") {
        chenx::dataset::median_impute(X);
    } else {
        throw std::invalid_argument("method must be `mean` or `median`");
    }
}

void encode(py::array_t<double>& X_arr, std::string_view method) {
    Mat<double> X = carma::arr_to_mat(X_arr, false);
    if (method == "dom") {
        chenx::dataset::dominance(X);
    } else if (method == "add" || method == "hybrid") {
    } else {
        throw std::invalid_argument("method must be `add`, `dom` or `hybrid`");
    }
}

py::array_t<double> hybird_value(const py::array_t<double>& genotype_arr, const py::array_t<double>& phenotype_arr) {
    Mat<double> genotype = carma::arr_to_mat_view(genotype_arr);
    Col<double> phenotype = carma::arr_to_col_view(phenotype_arr);
    return carma::mat_to_arr(chenx::dataset::hybird_value(genotype, phenotype));
}

void hybird(py::array_t<double>& genotype_arr, const py::array_t<double>& hybird_value_arr) {
    Mat<double> genotype = carma::arr_to_mat(genotype_arr);
    Mat<double> hybird_value = carma::arr_to_mat_view(hybird_value_arr);
    chenx::dataset::hybird(genotype, hybird_value);
}
}  // namespace binds

void bind_dataset(py::module& m) {
    m.def("_grm", &binds::grm, "Genetic relationship matrix", py::arg("X"), py::arg("encode_method") = "add");
    m.def("_hybrid_value", &binds::hybird_value, "calculate hybird value", py::arg("genotype"), py::arg("phenotype"));
    m.def("_hybrid", &binds::hybird, "Replace D value", py::arg("genotype"), py::arg("hybird_value"));
    m.def("_impute", &binds::impute, "Impute missing value", py::arg("X"), py::arg("method") = "mean");
    m.def("_encode", &binds::encode, "Encode missing value", py::arg("X"), py::arg("method") = "add");
}

PYBIND11_MODULE(_core, m) {
    bind_dataset(m);
}
