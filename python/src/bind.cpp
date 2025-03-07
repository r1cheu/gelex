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
#include <string_view>
#include <vector>

#include "array_caster.h"
#include "chenx/data/bed_reader.h"
#include "chenx/data/cross_grm.h"
#include "chenx/data/grm.h"
#include "chenx/estimator.h"
#include "chenx/model/linear_mixed_model.h"
#include "chenx/predictor.h"

namespace bind
{
namespace nb = nanobind;
using nb::literals::operator""_a;

NB_MODULE(_chenx, m)
{
    nb::class_<chenx::LinearMixedModelParams>(m, "_LinearMixedModelParams")
        .def(
            "__init__",
            [](chenx::LinearMixedModelParams* self,
               arr1d beta,
               arr1d sigma,
               arr1d proj_y,
               std::vector<std::string> dropped_individuals)
            {
                new (self) chenx::LinearMixedModelParams{
                    ToArma(beta),
                    ToArma(sigma),
                    ToArma(proj_y),
                    std::move(dropped_individuals)};
            },
            "beta"_a.noconvert(),
            "sigma"_a.noconvert(),
            "proj_y"_a.noconvert(),
            "dropped_individuals"_a.noconvert())
        .def(
            "__init__",
            [](chenx::LinearMixedModelParams* self,
               chenx::LinearMixedModel& model,
               std::vector<std::string> dropped_individuals)
            {
                new (self) chenx::LinearMixedModelParams{
                    model, std::move(dropped_individuals)};
            },
            "model"_a,
            "dropped_individuals"_a.noconvert())
        .def_prop_ro(
            "beta",
            [](const chenx::LinearMixedModelParams& self)
            { return ToPy(self.beta()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "sigma",
            [](const chenx::LinearMixedModelParams& self)
            { return ToPy(self.sigma()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "proj_y",
            [](const chenx::LinearMixedModelParams& self)
            { return ToPy(self.proj_y()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "dropped_individuals",
            [](const chenx::LinearMixedModelParams& self)
            { return self.dropped_individuals(); },
            nb::rv_policy::reference_internal);

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
            "num_fixed_effects",
            &chenx::LinearMixedModel::num_fixed_effects,
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "num_random_effects",
            &chenx::LinearMixedModel::num_random_effects,
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "num_individuals",
            &chenx::LinearMixedModel::num_individuals,
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "random_effect_names",
            &chenx::LinearMixedModel::random_effect_names,
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "_U",
            [](chenx::LinearMixedModel& self) { return ToPy(self.U()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "_proj_y",
            [](chenx::LinearMixedModel& self) { return ToPy(self.proj_y()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "beta",
            [](chenx::LinearMixedModel& self) { return ToPy(self.beta()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "sigma",
            [](chenx::LinearMixedModel& self) { return ToPy(self.sigma()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "y",
            [](chenx::LinearMixedModel& self) { return ToPy(self.y()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "X",
            [](chenx::LinearMixedModel& self) { return ToPy(self.X()); },
            nb::rv_policy::reference_internal)
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

    nb::class_<chenx::Predictor>(m, "_Predictor")
        .def(
            nb::init<std::string_view, chenx::LinearMixedModelParams&&>(),
            "train_bed"_a,
            "params"_a)
        .def(
            "set_cross_grm",
            [](chenx::Predictor& self,
               std::string_view method,
               arr1d center,
               double scale_factor,
               uint64_t chuck_size)
            {
                self.set_cross_grm(
                    method, ToRowVec(center), scale_factor, chuck_size);
            })
        .def(
            "_compute_u",
            [](chenx::Predictor& self, std::string_view test_bed)
            { return ToPy(self.ComputeU(test_bed)); })
        .def(
            "_compute_covariates",
            [](chenx::Predictor& self, arr2d covariates)
            { return ToPy(self.ComputeFixedEffects(ToArma(covariates))); })

        .def_prop_ro("test_individuals", &chenx::Predictor::test_individuals);

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
            "compute",
            [](chenx::AddGrm& self) { return ToPy(self.Compute()); },
            nb::rv_policy::move)
        .def_prop_ro(
            "individuals",
            [](chenx::AddGrm& self) { return self.bed().individuals(); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "center",
            [](chenx::AddGrm& self) { return ToPy(self.center()); },
            nb::rv_policy::reference_internal)
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
            "compute",
            [](chenx::DomGrm& self) { return ToPy(self.Compute()); },
            nb::rv_policy::move)
        .def_prop_ro(
            "individuals",
            [](chenx::DomGrm& self) { return self.bed().individuals(); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "center",
            [](chenx::DomGrm& self) { return ToPy(self.center()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "scale_factor",
            [](chenx::DomGrm& self) { return self.scale_factor(); });

    nb::class_<chenx::AddCrossGrm>(m, "add_cross_grm")
        .def(
            "__init__",
            [](chenx::AddCrossGrm* self,
               std::string_view train_bed_file,
               arr1d center,
               double scale_factor,
               uint64_t chunk_size,
               const std::vector<std::string>& exclude_individuals)
            {
                new (self) chenx::AddCrossGrm{
                    train_bed_file,
                    ToRowVec(center),
                    scale_factor,
                    chunk_size,
                    exclude_individuals};
            })
        .def(
            "compute",
            [](chenx::AddCrossGrm& self, std::string_view test_bed)
            { return ToPy(self.Compute(test_bed)); },
            nb::rv_policy::move);
}

}  // namespace bind
