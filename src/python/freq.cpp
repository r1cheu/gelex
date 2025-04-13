#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include <fmt/color.h>
#include <fmt/ranges.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>
#include <armadillo>

#include "gelex/data/bed_reader.h"
#include "gelex/data/grm.h"
#include "gelex/estimator/estimator.h"
#include "gelex/model/gblup.h"
#include "gelex/predictor.h"
#include "gelex/python/dense_caster.h"

namespace bind
{
namespace nb = nanobind;
using nb::literals::operator""_a;

void gblup_params(nb::module_& m)
{
    nb::class_<gelex::GBLUPParams>(m, "GBLUPParams")
        .def(
            "__init__",
            [](gelex::GBLUPParams* self,
               arr1d beta,
               arr1d sigma,
               arr1d proj_y,
               std::vector<std::string> dropped_individuals)
            {
                new (self) gelex::GBLUPParams{
                    to_arma(beta),
                    to_arma(sigma),
                    to_arma(proj_y),
                    std::move(dropped_individuals)};
            },
            "beta"_a.noconvert(),
            "sigma"_a.noconvert(),
            "proj_y"_a.noconvert(),
            "dropped_individuals"_a.noconvert(),
            nb::keep_alive<1, 2>(),
            nb::keep_alive<1, 3>(),
            nb::keep_alive<1, 4>())
        .def_prop_ro(
            "beta",
            [](const gelex::GBLUPParams& self) { return to_py(self.beta()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "sigma",
            [](const gelex::GBLUPParams& self) { return to_py(self.sigma()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "proj_y",
            [](const gelex::GBLUPParams& self) { return to_py(self.proj_y()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "dropped_individuals",
            [](const gelex::GBLUPParams& self)
            { return self.dropped_individuals(); },
            nb::rv_policy::reference_internal);
}

void gblup(nb::module_& m)
{
    nb::class_<gelex::GBLUP>(m, "_GBLUP")
        .def(
            "__init__",
            [](gelex::GBLUP* self, arr1d phenotype, arr2d design_mat_beta)
            {
                new (self)
                    gelex::GBLUP{to_arma(phenotype), to_arma(design_mat_beta)};
            },
            "phenotype"_a.noconvert(),
            "design_mat_beta"_a.noconvert(),
            nb::keep_alive<1, 2>(),
            nb::keep_alive<1, 3>())
        .def_prop_ro("n_common_effects", &gelex::GBLUP::n_common_effects)
        .def_prop_ro("n_group_effects", &gelex::GBLUP::n_genetic_effects)
        .def_prop_ro("n_random_effects", &gelex::GBLUP::n_random_effects)
        .def_prop_ro("n_individuals", &gelex::GBLUP::n_individuals)
        .def_prop_ro(
            "sigma_names",
            &gelex::GBLUP::sigma_names,
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "beta",
            [](gelex::GBLUP& self) { return to_py(self.beta()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "sigma",
            [](gelex::GBLUP& self) { return to_py(self.sigma()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "phenotype",
            [](gelex::GBLUP& self) { return to_py(self.phenotype()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "design_mat_beta",
            [](gelex::GBLUP& self) { return to_py(self.design_mat_beta()); },
            nb::rv_policy::reference_internal)
        .def(
            "add_genetic_effect",
            [](gelex::GBLUP& self, std::string name, arr2d& covar_mat)
            { self.add_genetic_effect(std::move(name), to_arma(covar_mat)); },
            nb::keep_alive<1, 3>())
        .def("reset", &gelex::GBLUP::reset, "reset the model")
        .def(
            "__repr__",
            [](const gelex::GBLUP& self)
            {
                return fmt::format(
                    "<GBLUP object at {:p}: {:d} Individuals, {:d} Common"
                    "effects, {:d} Group effects, {:d} Genetic effects, "
                    "Random effects: {}>",
                    static_cast<const void*>(&self),
                    self.n_individuals(),
                    self.n_common_effects(),
                    self.n_group_effects(),
                    self.n_genetic_effects(),
                    fmt::join(self.sigma_names(), ", "));
            })
        .def(
            "__str__",
            [](const gelex::GBLUP& self)
            {
                std::string info = fmt::format(
                    "┌─ GBLUP Model ─────────────────────────────────\n"
                    "│ Individuals:    {:6d}\n"
                    "│ Common Effects:    {:6d}\n"
                    "│ Group Effects:    {:6d}\n"
                    "│ Genetic Effects:  {:6d}\n"
                    "│ Sigma Names: {}\n"
                    "└───────────────────────────────────────────────",
                    self.n_individuals(),
                    self.n_common_effects(),
                    self.n_group_effects(),
                    self.n_genetic_effects(),
                    fmt::join(self.sigma_names(), ", "));

                return info;
            });
}

void predictor(nb::module_& m)
{
    nb::class_<gelex::Predictor>(m, "_Predictor")
        .def(
            nb::init<std::string_view, gelex::GBLUPParams>(),
            "train_bed"_a,
            "params"_a,
            nanobind::keep_alive<1, 3>())
        .def(
            "set_cross_grm",
            [](gelex::Predictor& self,
               std::string_view method,
               arr1d center,
               double scale_factor,
               uint64_t chuck_size)
            {
                self.set_cross_grm(
                    method, ToRowVec(center), scale_factor, chuck_size);
            },
            nb::keep_alive<1, 3>())
        .def(
            "_compute_random_effects",
            [](gelex::Predictor& self, std::string_view test_bed)
            { return to_py(self.compute_group_effects(test_bed)); })
        .def(
            "_compute_fixed_effects",
            [](gelex::Predictor& self, arr2d covariates)
            { return to_py(self.compute_common_effects(to_arma(covariates))); })

        .def_prop_ro("test_individuals", &gelex::Predictor::test_individuals);
}

void estimator(nb::module_& m)
{
    nb::class_<gelex::Estimator>(m, "Estimator")
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
            &gelex::Estimator::fit,
            "model"_a,
            "em_init"_a = true,
            "verbose"_a = true,
            "fit the model\n\n"
            "Parameters\n"
            "----------\n"
            "model : GBLUP\n"
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
            &gelex::Estimator::set_optimizer,
            "optimizer"_a = "NR",
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
}

void add_grm(nb::module_& m)
{
    nb::class_<gelex::AddGrm>(m, "add_grm")

        .def(
            nb::init<std::string_view, uint64_t, std::vector<std::string>>(),
            "bed_file"_a,
            "chunk_size"_a = gelex::DEFAULT_CHUNK_SIZE,
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
            [](gelex::AddGrm& self) { return to_py(self.compute()); },
            nb::rv_policy::move)
        .def_prop_ro(
            "individuals",
            [](gelex::AddGrm& self) { return self.bed().individuals(); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "center",
            [](gelex::AddGrm& self) { return to_py(self.center()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "scale_factor",
            [](gelex::AddGrm& self) { return self.scale_factor(); });
}

void dom_grm(nb::module_& m)
{
    nb::class_<gelex::DomGrm>(m, "dom_grm")
        .def(
            nb::init<std::string_view, uint64_t, std::vector<std::string>>(),
            "bed_file"_a,
            "chunk_size"_a = gelex::DEFAULT_CHUNK_SIZE,
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
            [](gelex::DomGrm& self) { return to_py(self.compute()); },
            nb::rv_policy::move)
        .def_prop_ro(
            "individuals",
            [](gelex::DomGrm& self) { return self.bed().individuals(); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "center",
            [](gelex::DomGrm& self) { return to_py(self.center()); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "scale_factor",
            [](gelex::DomGrm& self) { return self.scale_factor(); });
}

}  // namespace bind
