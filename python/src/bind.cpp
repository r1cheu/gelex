#include <fmt/color.h>
#include <fmt/ranges.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>
#include <armadillo>
#include <cstdint>
#include <string>
#include <vector>

#include "array_caster.h"
#include "chenx/data/bed_reader.h"
#include "chenx/data/grm.h"
#include "chenx/estimator.h"
#include "chenx/model/linear_mixed_model.h"

namespace bind
{
namespace nb = nanobind;
using nb::literals::operator""_a;

NB_MODULE(_chenx, m)
{
    nb::class_<chenx::LinearMixedModelParams>(m, "LinearMixedModelParams")
        .def(
            "__init__",
            [](chenx::LinearMixedModelParams* self,
               arr1d beta,
               arr1d sigma,
               std::vector<std::string> individuals,
               std::vector<std::string> dropped_individuals)
            {
                new (self) chenx::LinearMixedModelParams{
                    ToArma(beta),
                    ToArma(sigma),
                    std::move(individuals),
                    std::move(dropped_individuals)};
            },
            "beta"_a.noconvert(),
            "sigma"_a.noconvert(),
            "individuals"_a.noconvert(),
            "dropped_individuals"_a.noconvert())
        .def_prop_ro(
            "beta",
            [](chenx::LinearMixedModelParams& self) { return ToPy(self.beta); })
        .def_prop_ro(
            "sigma",
            [](chenx::LinearMixedModelParams& self)
            { return ToPy(self.sigma); })
        .def_prop_ro(
            "individuals",
            [](chenx::LinearMixedModelParams& self)
            { return self.individuals; })
        .def_prop_ro(
            "dropped_individuals",
            [](chenx::LinearMixedModelParams& self)
            { return self.dropped_individuals; });
    nb::class_<chenx::LinearMixedModel>(m, "_LinearMixedModel")
        .def(
            "__init__",
            [](chenx::LinearMixedModel* self,
               arr2d y,
               arr2d X,
               arr3d covar_mat,
               std::vector<std::string> names)
            {
                new (self) chenx::LinearMixedModel{
                    ToArma(y), ToArma(X), ToArma(covar_mat), std::move(names)};
            },
            "y"_a.noconvert(),
            "X"_a.noconvert(),
            "covar_mat"_a.noconvert(),
            "names"_a.noconvert())
        .def_prop_ro(
            "num_fixed_effects", &chenx::LinearMixedModel::num_fixed_effects)
        .def_prop_ro(
            "num_random_effects", &chenx::LinearMixedModel::num_random_effects)
        .def_prop_ro(
            "num_individuals", &chenx::LinearMixedModel::num_individuals)
        .def_prop_ro(
            "random_effect_names",
            &chenx::LinearMixedModel::random_effect_names)
        .def_prop_ro(
            "_U", [](chenx::LinearMixedModel& self) { return ToPy(self.U()); })
        .def_prop_ro(
            "beta",
            [](chenx::LinearMixedModel& self) { return ToPy(self.beta()); })
        .def_prop_ro(
            "sigma",
            [](chenx::LinearMixedModel& self) { return ToPy(self.sigma()); })
        .def("reset", &chenx::LinearMixedModel::Reset, "reset the model")
        .def(
            "__repr__",
            [](chenx::LinearMixedModel& self)
            {
                return fmt::format(
                    "Linear Mixed Model\n{:d} Individuals, {:d} Fixed effect, "
                    "Random Effect: [{}]",
                    self.num_individuals(),
                    self.num_fixed_effects(),
                    fmt::join(self.random_effect_names(), ", "));
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
            "verbose"_a = true,
            "Fit the model\n\n"
            "Parameters\n"
            "----------\n"
            "model : LinearMixedModel\n"
            "    The linear mixed model to fit\n"
            "em_init : bool, optional\n"
            "    Whether to use EM algorithm for initialization (default: "
            "True)\n\n"
            "verbose : bool, optional\n"
            "    Whether to print the optimization process (default: True)\n\n"
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
            nb::init<std::string_view, uint64_t, std::vector<std::string>>(),
            "bed_file"_a,
            "chunk_size"_a = chenx::DEFAULT_CHUNK_SIZE,
            "exclude_individuals"_a = std::vector<std::string>{},
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
            [](chenx::AddGrm& self) { return self.bed().individuals(); })
        .def_prop_ro(
            "center", [](chenx::AddGrm& self) { return ToPy(self.center()); })
        .def_prop_ro(
            "scale_factor",
            [](chenx::AddGrm& self) { return self.scale_factor(); });

    nb::class_<chenx::DomGrm>(m, "dom_grm")
        .def(
            nb::init<std::string_view, uint64_t, std::vector<std::string>>(),
            "bed_file"_a,
            "chunk_size"_a = chenx::DEFAULT_CHUNK_SIZE,
            "exclude_individuals"_a = std::vector<std::string>{},
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
            "center", [](chenx::DomGrm& self) { return ToPy(self.center()); })
        .def_prop_ro(
            "scale_factor",
            [](chenx::DomGrm& self) { return self.scale_factor(); });
}

}  // namespace bind
