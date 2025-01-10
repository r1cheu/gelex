#include <fmt/color.h>
#include <fmt/ranges.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>
#include <armadillo>

#include "array_caster.h"
#include "chenx/data/grm.h"
#include "chenx/estimator.h"
#include "chenx/model/linear_mixed_model.h"

namespace bind
{
NB_MODULE(_chenx, m)
{
    nb::class_<chenx::LinearMixedModel>(m, "LinearMixedModel")
        .def(
            "__init__",
            [](chenx::LinearMixedModel* self,
               const dvec& y,
               const dmat& X,
               const dcube& covar_mat,
               std::vector<std::string> names)
            {
                new (self) chenx::LinearMixedModel(
                    ToArma(y), ToArma(X), ToArma(covar_mat), std::move(names));
            },
            nb::arg("y").noconvert(),
            nb::arg("X").noconvert(),
            nb::arg("covar_mat").noconvert(),
            nb::arg("names"))
        .def(
            "__repr__",
            [](chenx::LinearMixedModel& self)
            {
                return fmt::format(
                    "Linear Mixed Model\n{:d} Samples, {:d} Fixed effect, "
                    "Random Effect: [{}]",
                    self.y().n_rows,
                    self.X().n_cols,
                    fmt::join(self.rand_names(), ", "));
            });
    nb::class_<chenx::Estimator>(m, "Estimator")
        .def(nb::init<std::string_view, size_t, double>())
        .def("fit", &chenx::Estimator::Fit)
        .def("set_optimizer", &chenx::Estimator::set_optimizer);
    m.def(
        "add_grm",
        [](dmat& genotype)
        {
            arma::dmat genotype_ = ToArma(genotype);
            return ToPy(chenx::AddGrm(genotype_));
        },
        nb::arg("genotype").noconvert(),
        nb::rv_policy::reference);
    m.def(
        "add_grm_chunk",
        [](dmat& genotype, dmat& grm)
        {
            arma::dmat genotype_ = ToArma(genotype);
            arma::dmat grm_ = ToArma(grm);
            chenx::AddGrmChunk(genotype_, grm_);
        },
        nb::arg("genotype").noconvert(),
        nb::arg("grm").noconvert());
    m.def(
        "dom_grm",
        [](dmat& genotype)
        {
            arma::dmat genotype_ = ToArma(genotype);
            return ToPy(chenx::DomGrm(genotype_));
        },
        nb::arg("genotype").noconvert(),
        nb::rv_policy::reference);
    m.def(
        "dom_grm_chunk",
        [](dmat& genotype, dmat& grm)
        {
            arma::dmat genotype_ = ToArma(genotype);
            arma::dmat grm_ = ToArma(grm);
            chenx::DomGrmChunk(genotype_, grm_);
        },
        nb::arg("genotype").noconvert(),
        nb::arg("grm").noconvert());
    m.def(
        "_scale_grm",
        [](dmat& grm)
        {
            arma::dmat grm_ = ToArma(grm);
            grm_ = grm_ / arma::trace(grm_) * static_cast<double>(grm_.n_rows);
        },
        nb::arg("grm").noconvert());
}

}  // namespace bind
