#include <fmt/format.h>
#include <fmt/ranges.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <armadillo>
#include <cstddef>

#include "dense.h"
#include "sparse.h"
#include "gelex/estimator/bayes/diagnostics.h"
#include "gelex/estimator/bayes/mcmc.h"
#include "gelex/estimator/bayes/result.h"
#include "gelex/estimator/bayes/samples.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/effects.h"

namespace bind
{
namespace nb = nanobind;
using nb::literals::operator""_a;
namespace gx = gelex;

using arma::dmat;
using arma::dvec;

void bayes_param(nb::module_& module)
{
    nb::class_<gx::SigmaPrior>(module, "SigmaParam")
        .def(nb::init<>())
        .def_rw("nu", &gx::SigmaPrior::nu)
        .def_rw("s2", &gx::SigmaPrior::s2);
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
    nb::class_<gx::BayesModel>(module, "BayesModel")
        .def(
            nb::init<std::string, dvec&&>(),
            "formula"_a,
            "phenotype"_a,
            nb::keep_alive<1, 3>())
        .def(
            "add_fixed_effect",
            &gx::BayesModel::add_fixed_effect,
            "names"_a,
            "levels"_a,
            "design_mat"_a,
            nb::keep_alive<1, 4>())
        .def(
            "add_random_effect",
            &gx::BayesModel::add_random_effect,
            "name"_a,
            "design_mat"_a,
            nb::keep_alive<1, 3>())
        .def(
            "add_genetic_effect",
            &gx::BayesModel::add_genetic_effect,
            "name"_a,
            "design_mat"_a,
            "type"_a,
            nb::keep_alive<1, 3>());
};

void mcmc_params(nb::module_& module)
{
    nb::class_<gx::MCMCParams>(module, "MCMCParams")
        .def(
            nb::init<size_t, size_t, size_t, size_t>(),
            "iter"_a,
            "n_burnin"_a,
            "n_thin"_a,
            "n_chains"_a)
        .def_rw("iter", &gx::MCMCParams::iter)
        .def_rw("n_burnin", &gx::MCMCParams::n_burnin)
        .def_rw("n_thin", &gx::MCMCParams::n_thin)
        .def_rw("n_chains", &gx::MCMCParams::n_chains);
}

void mcmc_storage(nb::module_& module)
{
    nb::class_<gx::MCMCSamples>(module, "MCMCSamples")
        .def_prop_ro("mu", &gx::MCMCSamples::mu)
        .def_prop_ro("fixed", &gx::MCMCSamples::fixed)
        .def_prop_ro("random", &gx::MCMCSamples::random)
        .def_prop_ro("genetic", &gx::MCMCSamples::genetic)
        .def_prop_ro("residual", &gx::MCMCSamples::residual);
}

void mcmc(nb::module_& module)
{
    nb::class_<gx::SampleGroup>(module, "SampleGroup")
        .def_rw(
            "coeffs",
            &gx::SampleGroup::coeffs,
            nb::rv_policy::reference_internal)
        .def_rw(
            "sigmas",
            &gx::SampleGroup::sigmas,
            nb::rv_policy::reference_internal);
    nb::class_<gx::MCMC>(module, "MCMC")
        .def(nb::init<gx::MCMCParams>(), "params"_a)
        .def("run", &gx::MCMC::run, "model"_a, "seed"_a = 42)
        .def("samples", &gx::MCMC::samples);
}

void mcmc_result(nb::module_& module)
{
    nb::class_<gx::PosteriorGroup>(module, "PosteriorGroup")
        .def_rw(
            "coeffs",
            &gx::PosteriorGroup::coeffs,
            nb::rv_policy::reference_internal)
        .def_rw(
            "sigmas",
            &gx::PosteriorGroup::sigmas,
            nb::rv_policy::reference_internal)
        .def_rw(
            "names",
            &gx::PosteriorGroup::names,
            nb::rv_policy::reference_internal);

    nb::class_<gx::PosteriorStats>(module, "PosteriorStats")
        .def("mean", &gx::PosteriorStats::mean, "index"_a = 0)
        .def("std", &gx::PosteriorStats::std, "index"_a = 0)
        .def("median", &gx::PosteriorStats::median, "index"_a = 0)
        .def("q5", &gx::PosteriorStats::q5, "index"_a = 0)
        .def("q95", &gx::PosteriorStats::q95, "index"_a = 0)
        .def("n_eff", &gx::PosteriorStats::n_eff, "index"_a = 0)
        .def("r_hat", &gx::PosteriorStats::r_hat, "index"_a = 0)

        .def_rw(
            "means",
            &gx::PosteriorStats::means,
            nb::rv_policy::reference_internal)
        .def_rw(
            "stds",
            &gx::PosteriorStats::stds,
            nb::rv_policy::reference_internal)
        .def_rw(
            "medians",
            &gx::PosteriorStats::medians,
            nb::rv_policy::reference_internal)
        .def_rw(
            "q5s", &gx::PosteriorStats::q5s, nb::rv_policy::reference_internal)
        .def_rw(
            "q95s",
            &gx::PosteriorStats::q95s,
            nb::rv_policy::reference_internal)
        .def_rw(
            "n_effs",
            &gx::PosteriorStats::n_effs,
            nb::rv_policy::reference_internal)
        .def_rw(
            "r_hats",
            &gx::PosteriorStats::r_hats,
            nb::rv_policy::reference_internal);

    nb::class_<gx::MCMCResult>(module, "MCMCResult")
        .def_rw("mu", &gx::MCMCResult::mu, nb::rv_policy::reference_internal)
        .def_rw(
            "fixed", &gx::MCMCResult::fixed, nb::rv_policy::reference_internal)
        .def_rw(
            "random",
            &gx::MCMCResult::random,
            nb::rv_policy::reference_internal)
        .def_rw(
            "genetic",
            &gx::MCMCResult::genetic,
            nb::rv_policy::reference_internal)
        .def_rw(
            "residual",
            &gx::MCMCResult::residual,
            nb::rv_policy::reference_internal);
}

void mcmc_diagnostics(nb::module_& m)
{
    m.def("gelman_rubin", &gelex::gelman_rubin, "samples"_a);
    m.def(
        "effective_sample_size",
        &gelex::effect_sample_size,
        "samples"_a,
        "bias"_a = true);
    m.def(
        "autocorrelation",
        &gelex::autocorrelation,
        "samples"_a,
        "bias"_a = true);
    m.def(
        "autocovariance", &gelex::autocovariance, "samples"_a, "bias"_a = true);
    m.def("split_gelman_rubin", &gelex::split_gelman_rubin, "samples"_a);
    m.def(
        "tran_cube",
        [](arma::dcube& x)
        {
            x.brief_print("before");
            x.slice(0) += 1;
            x.brief_print("after");
        });
}
}  // namespace bind
