#include <armadillo>
#include <cstdint>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>
#include <string>
#include <string_view>
#include <vector>

#include "array_caster.h"
#include "gelex/data/bed_reader.h"
#include "gelex/data/grm.h"
#include "gelex/estimator.h"
#include "gelex/model/linear_mixed_model.h"
#include "gelex/predictor.h"

namespace bind {
namespace nb = nanobind;
using nb::literals::operator""_a;

NB_MODULE(_gelex, m) {
  nb::class_<gelex::LinearMixedModelParams>(m, "_LinearMixedModelParams")
      .def(
          "__init__",
          [](gelex::LinearMixedModelParams *self, arr1d beta, arr1d sigma,
             arr1d proj_y, std::vector<std::string> dropped_individuals) {
            new (self) gelex::LinearMixedModelParams{
                ToArma(beta), ToArma(sigma), ToArma(proj_y),
                std::move(dropped_individuals)};
          },
          "beta"_a.noconvert(), "sigma"_a.noconvert(), "proj_y"_a.noconvert(),
          "dropped_individuals"_a.noconvert(), nb::keep_alive<1, 2>(),
          nb::keep_alive<1, 3>(), nb::keep_alive<1, 4>())
      .def(
          "__init__",
          [](gelex::LinearMixedModelParams *self,
             gelex::LinearMixedModel &model,
             std::vector<std::string> dropped_individuals) {
            new (self) gelex::LinearMixedModelParams{
                model, std::move(dropped_individuals)};
          },
          "model"_a, "dropped_individuals"_a.noconvert())
      .def_prop_ro(
          "beta",
          [](const gelex::LinearMixedModelParams &self) {
            return ToPy(self.beta());
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "sigma",
          [](const gelex::LinearMixedModelParams &self) {
            return ToPy(self.sigma());
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "proj_y",
          [](const gelex::LinearMixedModelParams &self) {
            return ToPy(self.proj_y());
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "dropped_individuals",
          [](const gelex::LinearMixedModelParams &self) {
            return self.dropped_individuals();
          },
          nb::rv_policy::reference_internal);

  nb::class_<gelex::LinearMixedModel>(m, "_LinearMixedModel")
      .def(
          "__init__",
          [](gelex::LinearMixedModel *self, arr2d y, arr2d X, arr3d covar_mat,
             std::vector<std::string> names) {
            new (self) gelex::LinearMixedModel{
                ToArma(y), ToArma(X), ToArma(covar_mat), std::move(names)};
          },
          "y"_a.noconvert(), "X"_a.noconvert(), "covar_mat"_a.noconvert(),
          "names"_a.noconvert(), nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>(),
          nb::keep_alive<1, 4>())
      .def_prop_ro("num_fixed_effects",
                   &gelex::LinearMixedModel::num_fixed_effects,
                   nb::rv_policy::reference_internal)
      .def_prop_ro("num_random_effects",
                   &gelex::LinearMixedModel::num_random_effects,
                   nb::rv_policy::reference_internal)
      .def_prop_ro("num_individuals", &gelex::LinearMixedModel::num_individuals,
                   nb::rv_policy::reference_internal)
      .def_prop_ro("random_effect_names",
                   &gelex::LinearMixedModel::random_effect_names,
                   nb::rv_policy::reference_internal)
      .def_prop_ro(
          "_U", [](gelex::LinearMixedModel &self) { return ToPy(self.U()); },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "_proj_y",
          [](gelex::LinearMixedModel &self) { return ToPy(self.proj_y()); },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "beta",
          [](gelex::LinearMixedModel &self) { return ToPy(self.beta()); },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "sigma",
          [](gelex::LinearMixedModel &self) { return ToPy(self.sigma()); },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "y", [](gelex::LinearMixedModel &self) { return ToPy(self.y()); },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "X", [](gelex::LinearMixedModel &self) { return ToPy(self.X()); },
          nb::rv_policy::reference_internal)
      .def("reset", &gelex::LinearMixedModel::Reset, "reset the model")
      .def("__repr__", [](gelex::LinearMixedModel &self) {
        return fmt::format(
            "Linear Mixed Model\n{:d} Individuals, {:d} Fixed effect, "
            "Random Effect: [{}]",
            self.num_individuals(), self.num_fixed_effects(),
            fmt::join(self.random_effect_names(), ", "));
      });

  nb::class_<gelex::Predictor>(m, "_Predictor")
      .def(nb::init<std::string_view, gelex::LinearMixedModelParams &&>(),
           "train_bed"_a, "params"_a, nanobind::keep_alive<1, 3>())
      .def(
          "set_cross_grm",
          [](gelex::Predictor &self, std::string_view method, arr1d center,
             double scale_factor, uint64_t chuck_size) {
            self.set_cross_grm(method, ToRowVec(center), scale_factor,
                               chuck_size);
          },
          nb::keep_alive<1, 3>())
      .def("_compute_random_effects",
           [](gelex::Predictor &self, std::string_view test_bed) {
             return ToPy(self.ComputeRandomEffects(test_bed));
           })
      .def("_compute_fixed_effects",
           [](gelex::Predictor &self, arr2d covariates) {
             return ToPy(self.ComputeFixedEffects(ToArma(covariates)));
           })

      .def_prop_ro("test_individuals", &gelex::Predictor::test_individuals);

  nb::class_<gelex::Estimator>(m, "Estimator")
      .def(nb::init<std::string_view, size_t, double>(), "optimizer"_a = "NR",
           "max_iter"_a = 20, "tol"_a = 1e-8,
           "Initialize the Estimator\n\n"
           "Parameters\n"
           "----------\n"
           "optimizer : str, optional\n"
           "    The optimization algorithm to use (default: 'NR')\n"
           "max_iter : int, optional\n"
           "    Maximum number of iterations (default: 20)\n"
           "tol : float, optional\n"
           "    Convergence tolerance (default: 1e-8)")
      .def("fit", &gelex::Estimator::Fit, "model"_a, "em_init"_a = true,
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
      .def("set_optimizer", &gelex::Estimator::set_optimizer,
           "optimizer"_a = "NR", "max_iter"_a = 20, "tol"_a = 1e-8,
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

  nb::class_<gelex::AddGrm>(m, "add_grm")

      .def(nb::init<std::string_view, uint64_t, std::vector<std::string>>(),
           "bed_file"_a, "chunk_size"_a = gelex::DEFAULT_CHUNK_SIZE,
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
          "compute", [](gelex::AddGrm &self) { return ToPy(self.Compute()); },
          nb::rv_policy::move)
      .def_prop_ro(
          "individuals",
          [](gelex::AddGrm &self) { return self.bed().individuals(); },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "center", [](gelex::AddGrm &self) { return ToPy(self.center()); },
          nb::rv_policy::reference_internal)
      .def_prop_ro("scale_factor",
                   [](gelex::AddGrm &self) { return self.scale_factor(); });

  nb::class_<gelex::DomGrm>(m, "dom_grm")
      .def(nb::init<std::string_view, uint64_t, std::vector<std::string>>(),
           "bed_file"_a, "chunk_size"_a = gelex::DEFAULT_CHUNK_SIZE,
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
          "compute", [](gelex::DomGrm &self) { return ToPy(self.Compute()); },
          nb::rv_policy::move)
      .def_prop_ro(
          "individuals",
          [](gelex::DomGrm &self) { return self.bed().individuals(); },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "center", [](gelex::DomGrm &self) { return ToPy(self.center()); },
          nb::rv_policy::reference_internal)
      .def_prop_ro("scale_factor",
                   [](gelex::DomGrm &self) { return self.scale_factor(); });
}

} // namespace bind
