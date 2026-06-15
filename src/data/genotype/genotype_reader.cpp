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

#include "gelex/data/genotype/genotype_reader.h"

#include <fmt/format.h>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype/genotype.h"
#include "gelex/data/genotype/method.h"
#include "gelex/data/locus_encoding.h"
#include "gelex/exception.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/geno_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/binary_writer.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::genotype
{

namespace
{

struct ChunkStats
{
    Eigen::VectorXd means;
    Eigen::VectorXd vars;
    Eigen::VectorXd A1freqs;
    std::vector<int64_t> mono_indices;
};

template <gelex::GeneticMode GT>
auto process_chunk(
    Eigen::Ref<Eigen::MatrixXd> chunk,
    Eigen::Index global_start,
    gelex::GenotypeMethod method,
    ChunkStats& stats) -> void
{
    const gelex::LociEncoding encoding{
        gelex::encode_inplace<double>(chunk, GT, method, 1e-12, global_start)};

    for (const auto& locus : encoding.loci)
    {
        stats.means[locus.marker_index] = locus.mean;
        stats.vars[locus.marker_index] = locus.var;
        stats.A1freqs[locus.marker_index]
            = locus.stats.has_nonmissing() ? locus.stats.A1freq() : 0.0;

        if (!locus.valid)
        {
            stats.mono_indices.push_back(
                static_cast<int64_t>(locus.marker_index));
        }
    }
}

auto load_mmapped(
    const std::filesystem::path& gbin_path,
    gelex::EffectType effect) -> MmappedStorage
{
    MmappedStorage mapped;
    mapped.reader
        = std::make_unique<gelex::io::BinaryReader>(gbin_path.string());

    auto geno_map
        = mapped.reader->to_map<double>(fmt::format("{}/genotype", effect));
    new (&mapped.view) MmappedStorage::MapType(
        geno_map.data(), geno_map.rows(), geno_map.cols());

    auto stats_mat
        = mapped.reader->to_mat<double>(fmt::format("{}/loci_stats", effect));
    mapped.mean = stats_mat.col(0);
    mapped.var = stats_mat.col(1);
    mapped.A1freq = stats_mat.col(2);

    if (const auto mono_path = fmt::format("{}/mono_indices", effect);
        mapped.reader->contains(mono_path))
    {
        auto mono_mat = mapped.reader->to_mat<int64_t>(mono_path);
        auto mono_col = mono_mat.col(0);
        mapped.mono_indices.resize(static_cast<size_t>(mono_col.size()));
        Eigen::Map<Eigen::Matrix<int64_t, Eigen::Dynamic, 1>>(
            mapped.mono_indices.data(),
            static_cast<Eigen::Index>(mapped.mono_indices.size())) = mono_col;
    }
    return mapped;
}

}  // namespace

GenotypeReader::GenotypeReader(
    const std::string& bfile_prefix,
    const dataframe::Index<std::string>& sample_index,
    gelex::GenoObserver observer)
    : bed_(gelex::open_bed(bfile_prefix, sample_index)),
      observer_(std::move(observer)),
      sample_size_(bed_.num_samples()),
      num_variants_(bed_.num_snps())
{
}

template <gelex::GeneticMode GT>
auto GenotypeReader::read(
    gelex::GenotypeMethod method,
    typename Sink::Variant sink,
    std::size_t chunk_size) -> Genotype
{
    return std::visit(
        [&](auto&& s) -> Genotype
        {
            using S = std::decay_t<decltype(s)>;
            if constexpr (std::is_same_v<S, Sink::InMemory>)
            {
                return read_in_memory<GT>(method, chunk_size);
            }
            else if constexpr (std::is_same_v<S, Sink::Mmap>)
            {
                return read_mmap<GT>(method, s.prefix, chunk_size);
            }
        },
        sink);
}

template <gelex::GeneticMode GT>
auto GenotypeReader::read_in_memory(
    gelex::GenotypeMethod method,
    std::size_t chunk_size) -> Genotype
{
    OwnedStorage owned;
    try
    {
        owned.data.resize(sample_size_, num_variants_);
    }
    catch (const std::bad_alloc&)
    {
        throw std::runtime_error(
            fmt::format(
                "Memory allocation failed for Genotype matrix ({} x {}). "
                "Requires approx {:.2f} GB RAM.",
                sample_size_,
                num_variants_,
                static_cast<double>(sample_size_)
                    * static_cast<double>(num_variants_) * sizeof(double)
                    / 1024.0 / 1024.0 / 1024.0));
    }

    ChunkStats stats;
    stats.means.resize(num_variants_);
    stats.vars.resize(num_variants_);
    stats.A1freqs.resize(num_variants_);
    stats.mono_indices.reserve(num_variants_ / 100);

    int64_t processed = 0;
    for (int64_t start = 0; start < num_variants_;)
    {
        const int64_t end
            = std::min(static_cast<int64_t>(start + chunk_size), num_variants_);
        auto target = owned.data.middleCols(start, end - start);
        bed_.read_into<double>(target, start);
        process_chunk<GT>(target, start, method, stats);
        processed += (end - start);
        gelex::notify(
            observer_,
            gelex::GenotypeProgressEvent{
                static_cast<size_t>(processed),
                static_cast<size_t>(num_variants_),
                false});
        start = end;
    }
    gelex::notify(
        observer_,
        gelex::GenotypeProgressEvent{
            static_cast<size_t>(num_variants_),
            static_cast<size_t>(num_variants_),
            true});

    std::ranges::sort(stats.mono_indices);
    owned.mono_indices = std::move(stats.mono_indices);
    owned.mean = std::move(stats.means);
    owned.var = std::move(stats.vars);
    owned.A1freq = std::move(stats.A1freqs);

    return Genotype(std::move(owned));
}

template <gelex::GeneticMode GT>
auto GenotypeReader::read_mmap(
    gelex::GenotypeMethod method,
    const std::filesystem::path& output_prefix,
    std::size_t chunk_size) -> Genotype
{
    constexpr auto effect = gelex::EffectType::from_genetic(GT);

    auto gbin_path = output_prefix;
    gbin_path += ".gbin";

    if (std::filesystem::exists(gbin_path))
    {
        auto logger = gelex::logging::get();
        logger->error("Output file already exists: [{}]", gbin_path.string());
        throw gelex::GelexException(
            fmt::format("{}: existing file", gbin_path.string()));
    }

    ChunkStats stats;
    stats.means.resize(num_variants_);
    stats.vars.resize(num_variants_);
    stats.A1freqs.resize(num_variants_);
    stats.mono_indices.reserve(num_variants_ / 100);

    {
        gelex::io::BinaryWriter writer(gbin_path.string());
        auto genotype_handle = writer.reserve<double>(
            fmt::format("{}/genotype", effect), sample_size_, num_variants_);

        int64_t processed = 0;
        for (int64_t start = 0; start < num_variants_;)
        {
            const int64_t end = std::min(
                static_cast<int64_t>(start + chunk_size), num_variants_);
            auto chunk = bed_.read<double>(start, end);
            process_chunk<GT>(chunk, start, method, stats);
            writer.write(genotype_handle, chunk);
            processed += (end - start);
            gelex::notify(
                observer_,
                gelex::GenotypeProgressEvent{
                    static_cast<size_t>(processed),
                    static_cast<size_t>(num_variants_),
                    false});
            start = end;
        }
        gelex::notify(
            observer_,
            gelex::GenotypeProgressEvent{
                static_cast<size_t>(num_variants_),
                static_cast<size_t>(num_variants_),
                true});

        std::ranges::sort(stats.mono_indices);

        auto stats_handle = writer.reserve<double>(
            fmt::format("{}/loci_stats", effect), num_variants_, 3);
        writer.write(stats_handle, stats.means);
        writer.write(stats_handle, stats.vars);
        writer.write(stats_handle, stats.A1freqs);

        if (!stats.mono_indices.empty())
        {
            auto mono_handle = writer.reserve<int64_t>(
                fmt::format("{}/mono_indices", effect),
                stats.mono_indices.size(),
                1);
            writer.write(mono_handle, stats.mono_indices);
        }
    }

    return Genotype(load_mmapped(gbin_path, effect));
}

auto GenotypeReader::read(
    const std::filesystem::path& gbin_path,
    gelex::GeneticMode mode) -> Genotype
{
    return Genotype(
        load_mmapped(gbin_path, gelex::EffectType::from_genetic(mode)));
}

template auto GenotypeReader::read<gelex::GeneticMode::A>(
    gelex::GenotypeMethod,
    typename Sink::Variant,
    std::size_t) -> Genotype;
template auto GenotypeReader::read<gelex::GeneticMode::D>(
    gelex::GenotypeMethod,
    typename Sink::Variant,
    std::size_t) -> Genotype;

}  // namespace gelex::genotype
