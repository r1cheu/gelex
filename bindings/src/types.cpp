#include <nanobind/nanobind.h>

#include "gelex/data/genotype_method.h"
#include "gelex/genetic_mode.h"

#include "gelex_py/register.h"

namespace nb = nanobind;

namespace gelex_py
{

void register_types(nb::module_& m)
{
    nb::enum_<gelex::GeneticMode>(m, "GeneticMode")
        .value("A", gelex::GeneticMode::A)
        .value("D", gelex::GeneticMode::D);

    nb::enum_<gelex::GenotypeMethod>(m, "GenotypeMethod")
        .value("StandardizeHWE", gelex::GenotypeMethod::StandardizeHWE)
        .value("CenterHWE", gelex::GenotypeMethod::CenterHWE)
        .value("Standardize", gelex::GenotypeMethod::Standardize)
        .value("Center", gelex::GenotypeMethod::Center)
        .value("OrthStandardizeHWE", gelex::GenotypeMethod::OrthStandardizeHWE)
        .value("OrthCenterHWE", gelex::GenotypeMethod::OrthCenterHWE)
        .value("OrthStandardize", gelex::GenotypeMethod::OrthStandardize)
        .value("OrthCenter", gelex::GenotypeMethod::OrthCenter)
        .value("NOIAStandardize", gelex::GenotypeMethod::NOIAStandardize)
        .value("NOIACenter", gelex::GenotypeMethod::NOIACenter);
}

}  // namespace gelex_py
