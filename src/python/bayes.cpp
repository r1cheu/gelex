#include <fmt/format.h>
#include <fmt/ranges.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <armadillo>

#include "dense.h"
#include "gelex/estimator/bayes/diagnostics.h"
#include "gelex/estimator/bayes/mcmc.h"
#include "gelex/estimator/bayes/result.h"
#include "gelex/estimator/bayes/samples.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"
#include "sparse.h"

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
            "_add_fixed_effect",
            &gx::BayesModel::add_fixed_effect,
            "names"_a,
            "levels"_a,
            "design_mat"_a,
            nb::keep_alive<1, 4>())
        .def(
            "_add_random_effect",
            &gx::BayesModel::add_random_effect,
            "name"_a,
            "design_mat"_a,
            nb::keep_alive<1, 3>())
        .def(
            "_add_genetic_effect",
            &gx::BayesModel::add_genetic_effect,
            "name"_a,
            "design_mat"_a,
            "type"_a,
            nb::keep_alive<1, 3>())
        .def("set_sigma_prior", &gx::BayesModel::set_sigma_prior)
        .def("set_sigma_prior_manual", &gx::BayesModel::set_sigma_prior_manul)
        .def("set_pi_prior", &gx::BayesModel::set_pi_prior)
        .def(
            "prior_summary",
            [](const gx::BayesModel& model)
            { nb::print(model.prior_summary().c_str()); })
        .def(
            "__repr__",
            [](const gx::BayesModel& model)
            {
                return fmt::format(
                    "<BayesModel(formula='{}') at {}>",
                    model.formula(),
                    static_cast<const void*>(&model));
            });
};

void mcmc_params(nb::module_& module)
{
    nb::class_<gx::MCMCParams>(module, "MCMCParams")
        .def(
            nb::init<size_t, size_t, size_t, size_t>(),
            "n_iters"_a = 5000,
            "n_burnin"_a = 3000,
            "n_thin"_a = 1,
            "n_chains"_a = 1,
            "Initialize MCMCParams with sampling parameters.\n\n"
            "    :param n_iters: Number of MCMC iterations (default: 5000).\n"
            "    :type n_iters: int\n"
            "    :param n_burnin: Number of burn-in iterations (default: "
            "3000).\n"
            "    :type n_burnin: int\n"
            "    :param n_thin: Thinning interval for samples (default: 1).\n"
            "    :type n_thin: int\n"
            "    :param n_chains: Number of independent chains (default: 1).\n"
            "    :type n_chains: int")
        .def_rw("n_iters", &gx::MCMCParams::n_iters)
        .def_rw("n_burnin", &gx::MCMCParams::n_burnin)
        .def_rw("n_thin", &gx::MCMCParams::n_thin)
        .def_rw("n_chains", &gx::MCMCParams::n_chains)
        .def(
            "__repr__",
            [](const gx::MCMCParams& params)
            {
                return fmt::format(
                    "<MCMCParams(n_iters={}, n_burnin={}, n_thin={}, "
                    "n_chains={}) at {}>",
                    params.n_iters,
                    params.n_burnin,
                    params.n_thin,
                    params.n_chains,
                    static_cast<const void*>(&params));
            })
        .def(
            "__str__",
            [](const gx::MCMCParams& params)
            {
                return fmt::format(
                    "MCMCParams: iters={}, burnin={}, thin={}, "
                    "chains={}",
                    params.n_iters,
                    params.n_burnin,
                    params.n_thin,
                    params.n_chains);
            });
}

void mcmc_storage(nb::module_& module)
{
    nb::class_<gx::MCMCSamples>(module, "MCMCSamples")
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
