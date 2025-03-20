#include <nanobind/nanobind.h>

namespace bind
{
namespace nb = nanobind;
using nb::literals::operator""_a;

void bayesa(nb::module_&);
void bayesrr(nb::module_&);
void bayesb(nb::module_&);
void bayesbpi(nb::module_&);
void bayesc(nb::module_&);
void bayescpi(nb::module_&);

void linear_mixed_model_params(nb::module_&);
void linear_mixed_model(nb::module_&);
void predictor(nb::module_&);
void estimator(nb::module_&);
void add_grm(nb::module_&);
void dom_grm(nb::module_&);

NB_MODULE(_core, m)  // NOLINT
{
    bayesa(m);
    bayesrr(m);
    bayesb(m);
    bayesbpi(m);
    bayesc(m);
    bayescpi(m);

    linear_mixed_model_params(m);
    linear_mixed_model(m);
    predictor(m);
    estimator(m);
    add_grm(m);
    dom_grm(m);
}
}  // namespace bind
