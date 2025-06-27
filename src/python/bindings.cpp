#include <nanobind/nanobind.h>

namespace bind
{
namespace nb = nanobind;
using nb::literals::operator""_a;
void sp_dense_dot(nb::module_& m);

void bayes_model(nb::module_&);
void bayes_param(nb::module_&);
void bayesalphabet(nb::module_&);

void mcmc_params(nb::module_&);
void mcmc(nb::module_&);
void mcmc_storage(nb::module_&);
void mcmc_result(nb::module_&);
void mcmc_diagnostics(nb::module_&);
void bayes_predictor(nb::module_&);

void gblup(nb::module_&);
void estimator(nb::module_&);
void bed_reader(nb::module_&);
void grm(nb::module_&);
void cross_grm(nb::module_&);

NB_MODULE(_core, m)  // NOLINT
{
    sp_dense_dot(m);

    bayesalphabet(m);
    bayes_param(m);
    bayes_model(m);
    bayes_predictor(m);

    mcmc_params(m);
    mcmc(m);
    mcmc_storage(m);
    mcmc_result(m);
    mcmc_diagnostics(m);

    grm(m);
    cross_grm(m);
    gblup(m);
    estimator(m);
    bed_reader(m);
}
}  // namespace bind
