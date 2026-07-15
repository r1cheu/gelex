/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "gelex/data/genotype_reader.h"

#include <Eigen/Core>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <memory>
#include <new>
#include <span>
#include <utility>
#include <vector>

#include "gelex/data/encode/encoder.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"
#include "gelex/data/genotype.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/snp_stats.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/geno_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/binary_writer.h"
#include "gelex/io/snpstats.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

namespace
{

struct ModePlan
{
    std::vector<gelex::GeneticMode> modes;
    std::vector<gelex::EncodingSpec> specs;
};

auto make_mode_plan(gelex::GeneticModeSet modes, gelex::GenotypeMethod method)
    -> ModePlan
{
    ModePlan plan;
    for (const gelex::GeneticMode mode : modes.each())
    {
        plan.modes.push_back(mode);
        plan.specs.push_back(gelex::encoding_spec_from_method(mode, method));
    }
    return plan;
}

auto make_stats(Eigen::Index num_variants) -> gelex::SnpStats
{
    gelex::SnpStats stats;
    stats.code.resize(3, num_variants);
    stats.A1freq.resize(num_variants);
    stats.valid_indices.reserve(num_variants);
    return stats;
}

// Fuses each variant straight from its packed form into every mode's chunk
// column: the genotype count is mode-independent, so it is tabulated once per
// variant and shared, and each mode's encoding then expands from it.
auto encode_chunk(
    const gelex::LocusEncoder& encoder,
    Eigen::Index start,
    std::span<const gelex::EncodingSpec> specs,
    std::span<Eigen::Ref<Eigen::MatrixXd>> targets,
    std::span<gelex::SnpStats> stats) -> void
{
    const Eigen::Index cols = targets.front().cols();
    for (Eigen::Index c = 0; c < cols; ++c)
    {
        const Eigen::Index marker = start + c;
        const gelex::LocusStats locus_stats{encoder.count(marker)};

        for (std::size_t k = 0; k < specs.size(); ++k)
        {
            const gelex::LocusEncoding encoding{
                encoder.encoding(marker, locus_stats, specs[k])};
            encoder.expand(marker, encoding, targets[k].col(c));

            stats[k].code.col(marker) = encoding.code.matrix();
            stats[k].A1freq[marker]
                = locus_stats.has_nonmissing() ? locus_stats.A1freq() : 0.0;
            if (encoding.valid)
            {
                stats[k].valid_indices.push_back(static_cast<int64_t>(marker));
            }
        }
    }
}

auto load_mmapped(
    const std::filesystem::path& geno_path,
    gelex::GeneticMode mode) -> MmappedStorage
{
    MmappedStorage mapped;
    mapped.reader = std::make_unique<gelex::BinaryReader>(geno_path.string());

    auto geno_map
        = mapped.reader->to_map<double>(fmt::format("{}/genotype", mode));
    new (&mapped.view) MmappedStorage::MapType(
        geno_map.data(), geno_map.rows(), geno_map.cols());

    mapped.stats = read_snp_stats(*mapped.reader, mode);
    return mapped;
}

}  // namespace

GenotypeReader::GenotypeReader(
    const gelex::Bed& bed,
    gelex::GenoObserver observer)
    : bed_(bed),
      observer_(std::move(observer)),
      sample_size_(bed_.num_samples()),
      num_variants_(bed_.num_snps())
{
}

auto GenotypeReader::read_in_memory(
    gelex::GeneticModeSet modes,
    gelex::GenotypeMethod method,
    std::size_t chunk_size) -> std::vector<ModeGenotype>
{
    const auto plan = make_mode_plan(modes, method);
    const std::size_t nm = plan.modes.size();

    std::vector<OwnedStorage> owned(nm);
    std::vector<SnpStats> stats;
    stats.reserve(nm);
    try
    {
        for (std::size_t k = 0; k < nm; ++k)
        {
            owned[k].data.resize(sample_size_, num_variants_);
            stats.push_back(make_stats(num_variants_));
        }
    }
    catch (const std::bad_alloc&)
    {
        const double gb = static_cast<double>(nm)
                          * static_cast<double>(sample_size_)
                          * static_cast<double>(num_variants_) * sizeof(double)
                          / 1024.0 / 1024.0 / 1024.0;
        throw gelex::GelexException(
            fmt::format(
                "Memory allocation failed for {} genotype matrix/matrices "
                "({} x {}). Requires approx {:.2f} GB RAM.",
                nm,
                sample_size_,
                num_variants_,
                gb));
    }

    const LocusEncoder encoder{bed_};
    int64_t processed = 0;
    for (int64_t start = 0; start < num_variants_;)
    {
        const int64_t end
            = std::min(start + static_cast<int64_t>(chunk_size), num_variants_);
        const int64_t cols = end - start;

        std::vector<Eigen::Ref<Eigen::MatrixXd>> targets;
        targets.reserve(nm);
        for (std::size_t k = 0; k < nm; ++k)
        {
            targets.emplace_back(owned[k].data.middleCols(start, cols));
        }
        encode_chunk(encoder, start, plan.specs, targets, stats);

        processed += cols;
        notify(
            observer_,
            GenotypeProgressEvent{
                static_cast<size_t>(processed),
                static_cast<size_t>(num_variants_),
                false});
        start = end;
    }
    notify(
        observer_,
        GenotypeProgressEvent{
            static_cast<size_t>(num_variants_),
            static_cast<size_t>(num_variants_),
            true});

    std::vector<ModeGenotype> result;
    result.reserve(nm);
    for (std::size_t k = 0; k < nm; ++k)
    {
        std::ranges::sort(stats[k].valid_indices);
        owned[k].stats = std::move(stats[k]);
        result.emplace_back(plan.modes[k], Genotype(std::move(owned[k])));
    }
    return result;
}

auto GenotypeReader::read_mmap(
    gelex::GeneticModeSet modes,
    gelex::GenotypeMethod method,
    const std::filesystem::path& output_prefix,
    std::size_t chunk_size) -> std::vector<ModeGenotype>
{
    const auto plan = make_mode_plan(modes, method);
    const std::size_t nm = plan.modes.size();

    std::vector<std::filesystem::path> paths(nm);
    std::vector<std::unique_ptr<BinaryWriter>> writers;
    std::vector<SectionHandle<double>> handles;
    std::vector<SnpStats> stats;
    writers.reserve(nm);
    handles.reserve(nm);
    stats.reserve(nm);
    for (std::size_t k = 0; k < nm; ++k)
    {
        auto geno_path = output_prefix;
        geno_path += fmt::format(".{}.geno", plan.modes[k]);
        if (std::filesystem::exists(geno_path))
        {
            throw gelex::GelexException(
                fmt::format(
                    "Output file already exists: [{}]", geno_path.string()));
        }
        paths[k] = geno_path;
        writers.push_back(std::make_unique<BinaryWriter>(geno_path.string()));
        handles.push_back(
            writers[k]->reserve<double>(
                fmt::format("{}/genotype", plan.modes[k]),
                sample_size_,
                num_variants_));
        stats.push_back(make_stats(num_variants_));
    }

    const LocusEncoder encoder{bed_};
    int64_t processed = 0;
    for (int64_t start = 0; start < num_variants_;)
    {
        const int64_t end
            = std::min(start + static_cast<int64_t>(chunk_size), num_variants_);
        const int64_t cols = end - start;

        std::vector<Eigen::MatrixXd> temps(
            nm, Eigen::MatrixXd(sample_size_, cols));
        std::vector<Eigen::Ref<Eigen::MatrixXd>> targets;
        targets.reserve(nm);
        for (auto& temp : temps)
        {
            targets.emplace_back(temp);
        }
        encode_chunk(encoder, start, plan.specs, targets, stats);
        for (std::size_t k = 0; k < nm; ++k)
        {
            writers[k]->write(handles[k], temps[k]);
        }

        processed += cols;
        notify(
            observer_,
            GenotypeProgressEvent{
                static_cast<size_t>(processed),
                static_cast<size_t>(num_variants_),
                false});
        start = end;
    }
    notify(
        observer_,
        GenotypeProgressEvent{
            static_cast<size_t>(num_variants_),
            static_cast<size_t>(num_variants_),
            true});

    for (std::size_t k = 0; k < nm; ++k)
    {
        std::ranges::sort(stats[k].valid_indices);
        write_snp_stats(*writers[k], plan.modes[k], stats[k]);
    }
    writers.clear();  // flush and close every file before mapping it back

    std::vector<ModeGenotype> result;
    result.reserve(nm);
    for (std::size_t k = 0; k < nm; ++k)
    {
        result.emplace_back(
            plan.modes[k], Genotype(load_mmapped(paths[k], plan.modes[k])));
    }
    return result;
}

auto GenotypeReader::read(
    const std::filesystem::path& geno_path,
    gelex::GeneticMode mode) -> Genotype
{
    return Genotype(load_mmapped(geno_path, mode));
}

}  // namespace gelex
