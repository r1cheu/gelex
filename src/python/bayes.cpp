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
#include "gelex/predictor/bayes/predictor.h"
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
    nb::class_<gx::ScaledInvChiSqParams>(module, "SigmaParam")
        .def(nb::init<>())
        .def_rw("nu", &gx::ScaledInvChiSqParams::nu)
        .def_rw("s2", &gx::ScaledInvChiSqParams::s2);
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
    nb::class_<gx::BayesModel>(module, "_BayesModel")
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
            "_mean",
            [](const gx::BayesModel& model) { return model.genetic()[0].mean; })
        .def(
            "_std",
            [](const gx::BayesModel& model)
            { return model.genetic()[0].stddev; })
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
    nb::class_<gx::MCMC>(module, "MCMC")
        .def(nb::init<gx::MCMCParams>(), "params"_a)
        .def("run", &gx::MCMC::run, "model"_a, "seed"_a = 42)
        .def("samples", &gx::MCMC::samples);
}

void mcmc_result(nb::module_& module)
{
    nb::class_<gx::PostieriorRandomSummary>(module, "PosteriorGroup")
        .def_rw(
            "coeffs",
            &gx::PostieriorRandomSummary::coeff,
            nb::rv_policy::reference_internal)
        .def_rw(
            "sigmas",
            &gx::PostieriorRandomSummary::sigma,
            nb::rv_policy::reference_internal);

    nb::class_<gx::PosteriorSummary>(module, "PosteriorSummary")
        .def_rw(
            "mean",
            &gx::PosteriorSummary::mean,
            nb::rv_policy::reference_internal)
        .def_rw(
            "std",
            &gx::PosteriorSummary::stddev,
            nb::rv_policy::reference_internal)
        .def_rw(
            "median",
            &gx::PosteriorSummary::median,
            nb::rv_policy::reference_internal)
        .def_rw(
            "hpdi_high",
            &gx::PosteriorSummary::hpdi_high,
            nb::rv_policy::reference_internal)
        .def_rw(
            "hpdi_low",
            &gx::PosteriorSummary::hpdi_low,
            nb::rv_policy::reference_internal)
        .def_rw(
            "ess",
            &gx::PosteriorSummary::ess,
            nb::rv_policy::reference_internal)
        .def_rw(
            "rhat",
            &gx::PosteriorSummary::rhat,
            nb::rv_policy::reference_internal);

    nb::class_<gx::MCMCResult>(module, "MCMCResult")
        .def_rw(
            "fixed", &gx::MCMCResult::fixed, nb::rv_policy::reference_internal)
        .def_rw(
            "snp_eff",
            &gx::MCMCResult::snp_eff,
            nb::rv_policy::reference_internal)
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
    m.def("hpdi", &gelex::hpdi, "samples"_a, "prob"_a = 0.90);
}

void bayes_predictor(nb::module_& m)
{
    nb::class_<gx::BayesPredictor>(m, "_BayesPredictor")
        .def(
            nb::init<const gx::BayesModel&, const gx::MCMCResult&>(),
            "model"_a,
            "result"_a)
        .def(
            "_predict",
            &gx::BayesPredictor::predict,
            "fixed_design"_a,
            "random_design"_a,
            "genetic_design"_a);
}
}  // namespace bind
