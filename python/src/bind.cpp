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
namespace nb = nanobind;
using nb::literals::operator""_a;
NB_MODULE(_chenx, m)
{
    nb::class_<chenx::LinearMixedModel>(m, "LinearMixedModel")
        .def(
            "__init__",
            [](chenx::LinearMixedModel* self,
               const dmat& y,
               const dmat& X,
               const dcube& covar_mat,
               std::vector<std::string> names)
            {
                new (self) chenx::LinearMixedModel(
                    ToArma(y), ToArma(X), ToArma(covar_mat), std::move(names));
            },
            "y"_a.noconvert(),
            "X"_a.noconvert(),
            "covar_mat"_a.noconvert(),
            "names"_a.noconvert())
        .def_prop_ro(
            "n_samples",
            [](chenx::LinearMixedModel& self) { return self.y().n_rows; })
        .def_prop_ro(
            "n_fixed_effect",
            [](chenx::LinearMixedModel& self) { return self.X().n_cols; })
        .def_prop_ro(
            "n_random_effect",
            [](chenx::LinearMixedModel& self)
            { return self.rand_names().size(); })
        .def_prop_ro(
            "beta",
            [](chenx::LinearMixedModel& self) { return ToPy(self.beta()); })
        .def_prop_ro(
            "sigma",
            [](chenx::LinearMixedModel& self) { return ToPy(self.sigma()); })
        .def_prop_ro(
            "_y", [](chenx::LinearMixedModel& self) { return ToPy(self.y()); })
        .def_prop_ro(
            "_X", [](chenx::LinearMixedModel& self) { return ToPy(self.X()); })
        .def_prop_ro(
            "_pdv",
            [](chenx::LinearMixedModel& self) { return ToPy(self.pdv()); })
        .def("reset", &chenx::LinearMixedModel::Reset, "reset the model")
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
        .def(
            nb::init<std::string_view, size_t, double>(),
            "optimizer"_a = "NR",
            "max_iter"_a = 20,
            "tol"_a = 1e-8,
            "Initialize the Estimator\n\n"
            "Parameters\n"
            "----------\n"
            "optimizer : str, optional\n"
            "    The optimization algorithm to use (default: 'NR')\n"
            "max_iter : int, optional\n"
            "    Maximum number of iterations (default: 20)\n"
            "tol : float, optional\n"
            "    Convergence tolerance (default: 1e-8)")
        .def(
            "fit",
            &chenx::Estimator::Fit,
            "model"_a,
            "em_init"_a = true,
            "Fit the model\n\n"
            "Parameters\n"
            "----------\n"
            "model : LinearMixedModel\n"
            "    The linear mixed model to fit\n"
            "em_init : bool, optional\n"
            "    Whether to use EM algorithm for initialization (default: "
            "True)\n\n"
            "Returns\n"
            "-------\n"
            "None")
        .def(
            "set_optimizer",
            &chenx::Estimator::set_optimizer,
            "optimizer"_a = "NR",
            "max_iter"_a = 20,
            "tol"_a = 1e-8,
            "reset optimizer\n\n"
            "Parameters\n"
            "----------\n"
            "optimizer : str, optional\n"
            "    The optimization algorithm to use (default: 'NR')\n"
            "max_iter : int, optional\n"
            "    Maximum number of iterations (default: 20)\n"
            "tol : float, optional\n"
            "    Convergence tolerance (default: 1e-8)\n\n"
            "Returns\n"
            "-------\n"
            "None");
    nb::class_<chenx::AddGrm>(m, "add_grm")
        .def(
            nb::init<const std::string&, uint64_t>(),
            "bed_file"_a,
            "chunk_size"_a = 10000,
            "Additive Genomic Relationship Matrix calculation.\n\n"
            "Parameters\n"
            "----------\n"
            "bed_file: str\n"
            "    The plink bed file path\n"
            "chunk_size: int, optional\n"
            "    Number of snps to processed per step (default: 10000)\n\n"
            "Returns\n"
            "-------\n"
            "np.ndarray\n")
        .def(
            "compute", [](chenx::AddGrm& self) { return ToPy(self.Compute()); })
        .def_prop_ro(
            "individuals",
            [](chenx::AddGrm& self) { return self.bed().individuals(); });

    nb::class_<chenx::DomGrm>(m, "dom_grm")
        .def(
            nb::init<const std::string&, uint64_t>(),
            "bed_file"_a,
            "chunk_size"_a = 10000,
            "Dominance Genomic Relationship Matrix calculation.\n\n"
            "Parameters\n"
            "----------\n"
            "bed_file: str\n"
            "    The plink bed file path\n"
            "chunk_size: int, optional\n"
            "    Number of snps to processed per step (default: 10000)\n\n"
            "Returns\n"
            "-------\n"
            "np.ndarray\n")
        .def(
            "compute", [](chenx::DomGrm& self) { return ToPy(self.Compute()); })
        .def_prop_ro(
            "individuals",
            [](chenx::DomGrm& self) { return self.bed().individuals(); })
        .def_prop_ro(
            "center", [](chenx::DomGrm& self) { return ToPy(self.center()); });
}

}  // namespace bind
