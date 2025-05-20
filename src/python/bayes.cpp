#include <fmt/format.h>
#include <fmt/ranges.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <armadillo>

#include "gelex/estimator/mcmc.h"
#include "gelex/estimator/mcmc_result.h"
#include "gelex/estimator/mcmc_storage.h"
#include "gelex/model/bayes.h"
#include "gelex/model/effects/bayes_effects.h"
#include "gelex/python/dense.h"
#include "gelex/python/sparse.h"

namespace bind
{
namespace nb = nanobind;
using nb::literals::operator""_a;
namespace gx = gelex;

using arma::dmat;
using arma::dvec;

void bayes_param(nb::module_& module)
{
    nb::class_<gx::bayes::SigmaParam>(module, "SigmaParam")
        .def(nb::init<>())
        .def_rw("nu", &gx::bayes::SigmaParam::nu)
        .def_rw("s2", &gx::bayes::SigmaParam::s2);
}

void bayesalphabet(nb::module_& module)
{
    nb::enum_<gx::BayesAlphabet>(module, "BayesAlphabet")
        .value("RR", gx::BayesAlphabet::RR)
        .value("A", gx::BayesAlphabet::A)
        .value("B", gx::BayesAlphabet::B)
        .value("Bpi", gx::BayesAlphabet::Bpi)
        .value("C", gx::BayesAlphabet::C)
        .value("Cpi", gx::BayesAlphabet::Cpi)
        .value("R", gx::BayesAlphabet::R)
        .export_values();
}

void bayes_model(nb::module_& module)
{
    nb::class_<gx::Bayes>(module, "Bayes")
        .def(
            nb::init<std::string, dvec&&>(),
            "formula"_a,
            "phenotype"_a,
            nb::keep_alive<1, 2>(),
            nb::keep_alive<1, 3>())
        .def(
            "add_fixed_effect",
            &gx::Bayes::add_fixed_effect,
            "names"_a,
            "levels"_a,
            "design_mat"_a,
            nb::keep_alive<1, 4>())
        .def(
            "add_random_effect",
            &gx::Bayes::add_random_effect,
            "name"_a,
            "design_mat"_a,
            nb::keep_alive<1, 3>())
        .def(
            "add_genetic_effect",
            &gx::Bayes::add_genetic_effect,
            "name"_a,
            "design_mat"_a,
            "type"_a,
            nb::keep_alive<1, 3>())
        .def(
            "__str__",
            [](const gx::Bayes& self)
            {
                return fmt::format(
                    "┌─ Bayes Model ─────────────────────────────────\n"
                    "│ Individuals:    {:6d}\n"
                    "│ Common Effects:    {}\n"
                    "│ Random Effects:    {}\n"
                    "│ Genetic Effects:  {}\n"
                    "└───────────────────────────────────────────────",
                    self.n_individuals(),
                    fmt::join(self.fixed().levels, ", "),
                    fmt::join(self.random().names(), ", "),
                    fmt::join(self.genetic().names(), ", "));
            });
}

void mcmc_params(nb::module_& module)
{
    nb::class_<gx::MCMCParams>(module, "MCMCParams")
        .def(
            nb::init<size_t, size_t, size_t, size_t>(),
            "iter"_a,
            "n_burnin"_a,
            "n_thin"_a,
            "seed"_a)
        .def_rw("iter", &gx::MCMCParams::iter)
        .def_rw("n_burnin", &gx::MCMCParams::n_burnin)
        .def_rw("n_thin", &gx::MCMCParams::n_thin)
        .def_rw("seed", &gx::MCMCParams::seed);
}

void mcmc_storage(nb::module_& module)
{
    nb::class_<gx::MCMCStorage>(module, "MCMCStorage")
        .def_prop_ro("mu_samples", &gx::MCMCStorage::mu_samples)
        .def_prop_ro("fixed_samples", &gx::MCMCStorage::fixed_samples)
        .def_prop_ro("random_samples", &gx::MCMCStorage::random_samples)
        .def_prop_ro("genetic_samples", &gx::MCMCStorage::genetic_samples)
        .def_prop_ro("residual_samples", &gx::MCMCStorage::residual_samples)
        .def_prop_ro(
            "random_sigma_samples", &gx::MCMCStorage::random_sigma_samples)
        .def_prop_ro(
            "genetic_sigma_samples", &gx::MCMCStorage::genetic_sigma_samples);
}

void mcmc(nb::module_& module)
{
    nb::class_<gx::MCMC>(module, "MCMC")
        .def(nb::init<gx::MCMCParams>(), "params"_a)
        .def("run", &gx::MCMC::run, "model"_a, "log_freq"_a = 100)
        .def("result", &gx::MCMC::result, nb::rv_policy::reference_internal)
        .def("storage", &gx::MCMC::storage, nb::rv_policy::reference_internal);
}

void mcmc_result(nb::module_& module)
{
    nb::class_<gx::ParameterResult>(module, "ParameterResult")
        .def_ro("mean", &gx::ParameterResult::mean)
        .def_ro("std", &gx::ParameterResult::std)
        .def_ro("median", &gx::ParameterResult::median)
        .def_ro("q5", &gx::ParameterResult::q5)
        .def_ro("q95", &gx::ParameterResult::q95)
        .def_ro("n_eff", &gx::ParameterResult::n_eff)
        .def_ro("r_hat", &gx::ParameterResult::r_hat);

    nb::class_<gx::EffectResult>(module, "EffectResult")
        .def("mean", &gx::EffectResult::mean, "index"_a = 0)
        .def("std", &gx::EffectResult::std, "index"_a = 0)
        .def("median", &gx::EffectResult::median, "index"_a = 0)
        .def("q5", &gx::EffectResult::q5, "index"_a = 0)
        .def("q95", &gx::EffectResult::q95, "index"_a = 0)
        .def("n_eff", &gx::EffectResult::n_eff, "index"_a = 0)
        .def("r_hat", &gx::EffectResult::r_hat, "index"_a = 0)
        .def("means", &gx::EffectResult::means)
        .def("stds", &gx::EffectResult::stds)
        .def("medians", &gx::EffectResult::medians)
        .def("q5s", &gx::EffectResult::q5s)
        .def("q95s", &gx::EffectResult::q95s)
        .def("n_effs", &gx::EffectResult::n_effs)
        .def("r_hats", &gx::EffectResult::r_hats);

    nb::class_<gx::MCMCResult>(module, "MCMCResult")
        .def_ro("mu", &gx::MCMCResult::mu)
        .def_ro("fixed", &gx::MCMCResult::fixed)
        .def_ro("random", &gx::MCMCResult::random)
        .def_ro("genetic", &gx::MCMCResult::genetic)
        .def_ro("residual", &gx::MCMCResult::residual)
        .def_ro("random_sigma", &gx::MCMCResult::random_sigma)
        .def_ro("genetic_sigma", &gx::MCMCResult::genetic_sigma)
        .def_ro("random_names", &gx::MCMCResult::random_names)
        .def_ro("genetic_names", &gx::MCMCResult::genetic_names);
}
}  // namespace bind
