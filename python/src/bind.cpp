#include <fmt/color.h>
#include <fmt/ranges.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>
#include <armadillo>

#include "arma.h"
#include "chenx/data/grm.h"
#include "chenx/estimator.h"
#include "chenx/model/linear_mixed_model.h"

namespace bind
{
namespace nb = nanobind;

class Test
{
   public:
    Test() = default;
    const arma::dmat& mat() const { return mat_; }

   private:
    arma::dmat mat_{3, 3, arma::fill::eye};
};

// note, the nb::ndarray is like a shared_ptr, so we choose to pass by value
// and add const to indicate that we don't want to modify the data

NB_MODULE(_chenx, m)
{
    nb::class_<Test>(m, "Test")
        .def(nb::init<>())
        .def_prop_ro(
            "mat",
            [](Test& t) { return ToPyView<dmat_view>(t.mat()); },
            nb::rv_policy::reference_internal);

    nb::class_<chenx::LinearMixedModel>(m, "LinearMixedModel")
        .def(
            "__init__",
            [](chenx::LinearMixedModel* t,
               dvec y,
               dmat_a X,
               dcube covar_mat,
               std::vector<std::string> names)
            {
                new (t) chenx::LinearMixedModel(
                    arma::dvec(y.data(), y.shape(0), false, true),
                    arma::dmat(X.data(), X.shape(0), X.shape(1), false, true),
                    arma::dcube(
                        covar_mat.data(),
                        covar_mat.shape(0),
                        covar_mat.shape(1),
                        covar_mat.shape(2),
                        false,
                        true),
                    std::move(names));
            },
            nb::arg("y").noconvert(),
            nb::arg("X").noconvert(),
            nb::arg("covar_mat").noconvert(),
            nb::arg("names"))
        .def(
            "__repr__",
            [](chenx::LinearMixedModel& t)
            {
                return fmt::format(
                    "Linear Mixed Model\n{:d} Samples, {:d} Fixed effect, "
                    "Random Effect: [{}]",
                    t.y().n_rows,
                    t.X().n_cols,
                    fmt::join(t.rand_names(), ", "));
            });
    nb::class_<chenx::Estimator>(m, "Estimator")
        .def(nb::init<std::string_view, size_t, double>())
        .def("fit", &chenx::Estimator::Fit)
        .def("set_optimizer", &chenx::Estimator::set_optimizer);
    m.def(
        "return_arma",
        []()
        {
            arma::dmat mat(3, 3, arma::fill::eye);
            return ToPy<dmat>(std::move(mat));
        });
    m.def(
        "compute_add_grm",
        [](dmat genotype)
        {
            arma::dmat genotype_mat(
                genotype.data(),
                genotype.shape(0),
                genotype.shape(1),
                false,
                true);

            return ToPy<dmat>(chenx::AdditiveGrm(genotype_mat));
        },
        nb::rv_policy::reference);
    m.def(
        "compute_dom_grm",
        [](dmat genotype)
        {
            arma::dmat genotype_mat(
                genotype.data(),
                genotype.shape(0),
                genotype.shape(1),
                false,
                true);
            return ToPy<dmat>(chenx::DomainanceGrm(genotype_mat));
        },
        nb::rv_policy::reference);
}

}  // namespace bind
