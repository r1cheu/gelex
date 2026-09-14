// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <filesystem>
#include <fmt/format.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/genotype/gebv.h"
#include "gelex/bayes/genotype/projection.h"
#include "gelex/data/bed.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/sample_id.h"
#include "gelex/data/snp_lut.h"
#include "gelex/data/snp_lut_io.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

#include "gelex_py/register.h"

namespace nb = nanobind;

namespace gelex_py
{

namespace
{

using gelex::bayes::CompactGenotype;
using gelex::bayes::GeneticDesign;
using gelex::bayes::GeneticProjection;

auto keys_of(const gelex::DataFrameIndex<std::string>& index)
    -> std::vector<std::string>
{
    return {index.keys().begin(), index.keys().end()};
}

using SampleId = std::pair<std::string, std::string>;

// The core keys samples as "FID<US>IID"; Python sees (fid, iid) pairs.
auto sample_ids_of(const gelex::Bed& bed) -> std::vector<SampleId>
{
    std::vector<SampleId> ids;
    ids.reserve(bed.sample_index().size());
    for (const auto& key : bed.sample_index().keys())
    {
        const auto [fid, iid] = gelex::split_sample_id(key);
        ids.emplace_back(fid, iid);
    }
    return ids;
}

auto sample_keys_of(const std::vector<SampleId>& ids)
    -> std::vector<std::string>
{
    std::vector<std::string> keys;
    keys.reserve(ids.size());
    for (const auto& [fid, iid] : ids)
    {
        if (fid.empty() || iid.empty())
        {
            throw gelex::GelexException(
                "gather: FID and IID must both be non-empty");
        }
        keys.push_back(
            fmt::format("{}{}{}", fid, gelex::sample_id_separator, iid));
    }
    return keys;
}

auto modes_of(const GeneticDesign& design) -> std::vector<gelex::GeneticMode>
{
    std::vector<gelex::GeneticMode> modes;
    for (const auto mode : design.each_mode())
    {
        modes.push_back(mode);
    }
    return modes;
}

auto design_from_luts(
    const gelex::Bed& bed,
    const gelex::ModeMap<gelex::SnpLutMatrix>& luts) -> GeneticDesign
{
    return gelex::bayes::make_genetic_design(bed, luts);
}

auto gebv_of(
    const GeneticProjection& projection,
    const Eigen::Ref<const Eigen::VectorXd>& coefficients) -> Eigen::VectorXd
{
    if (coefficients.size() != projection.cols())
    {
        throw gelex::GelexException(
            fmt::format(
                "gebv: {} coefficients for {} markers",
                coefficients.size(),
                projection.cols()));
    }
    const Eigen::Map<const Eigen::MatrixXd> as_matrix{
        coefficients.data(), coefficients.size(), 1};
    Eigen::VectorXd result(projection.rows());
    gelex::gebv_draw(projection, as_matrix, 0, result);
    return result;
}

auto register_bed(nb::module_& m) -> void
{
    nb::class_<gelex::Bed>(
        m,
        "Bed",
        "A PLINK1 .bed dataset with its .fam/.bim metadata. gather() narrows "
        "and reorders the samples; every design built afterwards follows that "
        "order.")
        .def_prop_ro("num_samples", &gelex::Bed::num_samples)
        .def_prop_ro("num_markers", &gelex::Bed::num_snps)
        .def_prop_ro(
            "sample_ids",
            &sample_ids_of,
            "(FID, IID) of every sample in the current order.")
        .def_prop_ro(
            "marker_ids",
            [](const gelex::Bed& bed) { return keys_of(bed.snp_index()); })
        .def(
            "gather",
            [](gelex::Bed& bed, const std::vector<SampleId>& samples)
            {
                bed.gather(
                    gelex::DataFrameIndex<std::string>{
                        sample_keys_of(samples)});
            },
            nb::arg("samples"),
            "Restrict the dataset to the given (FID, IID) pairs, in that "
            "order. Every pair must exist in the .fam file.");

    m.def(
        "open_bed",
        &gelex::open_bed,
        nb::arg("bfile"),
        "Open `bfile`.bed/.bim/.fam.");
}

auto register_projection(nb::module_& m) -> void
{
    nb::class_<GeneticProjection>(
        m,
        "GeneticProjection",
        "One genetic mode's encoded view of a design's genotypes: each marker "
        "is a 4-entry lookup table over the raw BED codes.")
        .def_prop_ro("num_samples", &GeneticProjection::rows)
        .def_prop_ro("num_markers", &GeneticProjection::cols)
        .def_prop_ro(
            "snp_luts",
            [](const GeneticProjection& projection) -> Eigen::MatrixXd
            { return projection.snp_luts().matrix(); },
            "(4, markers) lookup tables indexed by raw BED code "
            "(A1A1, missing, A1A2, A2A2).");

    m.def(
        "gebv",
        &gebv_of,
        nb::arg("projection"),
        nb::arg("coefficients"),
        "Genomic values of every sample under `projection` for one vector of "
        "marker coefficients (one posterior draw).");
}

auto register_design(nb::module_& m) -> void
{
    nb::class_<GeneticDesign>(
        m,
        "GeneticDesign",
        "The genotypes of one Bed's current samples together with a "
        "projection per genetic mode.")
        .def_prop_ro("num_samples", &GeneticDesign::rows)
        .def_prop_ro("num_markers", &GeneticDesign::cols)
        .def_prop_ro("modes", &modes_of)
        .def_prop_ro(
            "a1_frequency",
            &GeneticDesign::a1_frequency,
            nb::rv_policy::reference_internal)
        .def("contains", &GeneticDesign::contains, nb::arg("mode"))
        .def(
            "projection",
            &GeneticDesign::projection,
            nb::arg("mode"),
            nb::rv_policy::reference_internal);

    m.def(
        "make_genetic_design",
        &design_from_luts,
        nb::arg("bed"),
        nb::arg("luts"),
        "Build the design from the lookup tables of a trained model, as "
        "returned by load_snp_luts(). The Bed must carry the same markers in "
        "the same order and allele orientation as the training data.");
    m.def(
        "load_snp_luts",
        [](const std::string& path)
        { return gelex::load_snp_luts(std::filesystem::path{path}); },
        nb::arg("path"),
        "Read a .snplut file into {GeneticMode: (4, markers) array}.");
    m.def(
        "write_snp_luts",
        [](const std::string& path,
           const gelex::ModeMap<gelex::SnpLutMatrix>& luts)
        { gelex::write_snp_luts(std::filesystem::path{path}, luts); },
        nb::arg("path"),
        nb::arg("luts"));
}

}  // namespace

void register_genotype(nb::module_& m)
{
    register_bed(m);
    register_projection(m);
    register_design(m);
}

}  // namespace gelex_py
