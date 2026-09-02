#include "gelex/data/encode/types.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/vector.h>

#include "gelex/data/encode/stats.h"
#include "gelex/types/genetic_mode.h"

#include "gelex_py/register.h"

namespace nb = nanobind;

namespace gelex_py
{

void register_types(nb::module_& m)
{
    using namespace gelex;

    nb::enum_<GeneticMode>(m, "GeneticMode")
        .value("A", GeneticMode::A)
        .value("D", GeneticMode::D);

    nb::enum_<GenotypeMethod>(m, "GenotypeMethod")
        .value("StandardizeHWE", GenotypeMethod::StandardizeHWE)
        .value("CenterHWE", GenotypeMethod::CenterHWE)
        .value("Standardize", GenotypeMethod::Standardize)
        .value("Center", GenotypeMethod::Center)
        .value("OrthStandardizeHWE", GenotypeMethod::OrthStandardizeHWE)
        .value("OrthCenterHWE", GenotypeMethod::OrthCenterHWE)
        .value("OrthStandardize", GenotypeMethod::OrthStandardize)
        .value("OrthCenter", GenotypeMethod::OrthCenter)
        .value("NOIAStandardize", GenotypeMethod::NOIAStandardize)
        .value("NOIACenter", GenotypeMethod::NOIACenter);

    nb::enum_<DominanceCode>(m, "DominanceCode")
        .value("Het", DominanceCode::Het)
        .value("HWE", DominanceCode::HWE)
        .value("NOIA", DominanceCode::NOIA);

    nb::enum_<Normalization>(m, "Normalization")
        .value("None", Normalization::None)
        .value("Center", Normalization::Center)
        .value("CenterScale", Normalization::CenterScale);

    nb::enum_<MomentBasis>(m, "MomentBasis")
        .value("Empirical", MomentBasis::Empirical)
        .value("Theoretical", MomentBasis::Theoretical);

    nb::class_<LocusStats>(m, "LocusStats")
        .def_ro("nA2A2", &LocusStats::nA2A2)
        .def_ro("nA1A2", &LocusStats::nA1A2)
        .def_ro("nA1A1", &LocusStats::nA1A1)
        .def_ro("n_missing", &LocusStats::n_missing)
        .def("n_nonmissing", &LocusStats::n_nonmissing)
        .def("has_nonmissing", &LocusStats::has_nonmissing)
        .def("pA2A2", &LocusStats::pA2A2)
        .def("pA1A2", &LocusStats::pA1A2)
        .def("pA1A1", &LocusStats::pA1A1)
        .def("A1freq", &LocusStats::A1freq);

    nb::class_<EncodingSpec>(m, "EncodingSpec")
        .def_ro("effect", &EncodingSpec::effect)
        .def_ro("dominance_code", &EncodingSpec::dominance_code)
        .def_ro("normalization", &EncodingSpec::normalization)
        .def_ro("moment_basis", &EncodingSpec::moment_basis);

    nb::class_<LocusEncoding>(m, "LocusEncoding")
        .def_ro("column_index", &LocusEncoding::column_index)
        .def_ro("marker_index", &LocusEncoding::marker_index)
        .def_ro("stats", &LocusEncoding::stats)
        .def_ro("lut", &LocusEncoding::lut)
        .def_ro("mean", &LocusEncoding::mean)
        .def_ro("var", &LocusEncoding::var)
        .def_ro("sd", &LocusEncoding::sd)
        .def_ro("valid", &LocusEncoding::valid);
}

}  // namespace gelex_py
