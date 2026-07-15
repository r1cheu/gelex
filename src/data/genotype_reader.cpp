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
#include <span>
#include <utility>
#include <vector>

#include "gelex/data/encode/encoder.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"
#include "gelex/data/genotype.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/genotype_sink.h"
#include "gelex/data/snp_stats.h"
#include "gelex/infra/logging/geno_event.h"
#include "gelex/infra/logging/notify.h"
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

// Single mode-agnostic pass: count each variant once, expand into whatever
// columns the sink hands out, then let the sink persist the chunk. The sink is
// the only thing that differs between resident and mmap output.
auto run_encode(
    const gelex::Bed& bed,
    const gelex::GenoObserver& observer,
    const ModePlan& plan,
    gelex::GenotypeSink& sink,
    std::size_t chunk_size) -> gelex::ModeMap<gelex::Genotype>
{
    const int64_t num_variants = bed.num_snps();
    const std::size_t nm = plan.modes.size();

    std::vector<gelex::SnpStats> stats;
    stats.reserve(nm);
    for (std::size_t k = 0; k < nm; ++k)
    {
        stats.push_back(make_stats(num_variants));
    }

    const gelex::LocusEncoder encoder{bed};
    int64_t processed = 0;
    for (int64_t start = 0; start < num_variants;)
    {
        const int64_t end
            = std::min(start + static_cast<int64_t>(chunk_size), num_variants);
        const int64_t cols = end - start;

        auto targets = sink.chunk_targets(start, cols);
        encode_chunk(encoder, start, plan.specs, targets, stats);
        sink.commit_chunk();

        processed += cols;
        notify(
            observer,
            gelex::GenotypeProgressEvent{
                static_cast<size_t>(processed),
                static_cast<size_t>(num_variants),
                false});
        start = end;
    }
    notify(
        observer,
        gelex::GenotypeProgressEvent{
            static_cast<size_t>(num_variants),
            static_cast<size_t>(num_variants),
            true});

    for (auto& s : stats)
    {
        std::ranges::sort(s.valid_indices);
    }
    return sink.finalize(stats);
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
    std::size_t chunk_size) -> gelex::ModeMap<gelex::Genotype>
{
    const auto plan = make_mode_plan(modes, method);
    InMemorySink sink{plan.modes, sample_size_, num_variants_};
    return run_encode(bed_, observer_, plan, sink, chunk_size);
}

auto GenotypeReader::read_mmap(
    gelex::GeneticModeSet modes,
    gelex::GenotypeMethod method,
    const std::filesystem::path& output_prefix,
    std::size_t chunk_size) -> gelex::ModeMap<gelex::Genotype>
{
    const auto plan = make_mode_plan(modes, method);
    MmapSink sink{plan.modes, sample_size_, num_variants_, output_prefix};
    return run_encode(bed_, observer_, plan, sink, chunk_size);
}

auto GenotypeReader::read(
    const std::filesystem::path& geno_path,
    gelex::GeneticMode mode) -> Genotype
{
    return load_genotype(geno_path, mode);
}

}  // namespace gelex
