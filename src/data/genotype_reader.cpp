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

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/locus_encoding.h"
#include "gelex/exception.h"
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

auto process_chunk(
    Eigen::Ref<Eigen::MatrixXd> chunk,
    Eigen::Index global_start,
    gelex::GeneticMode mode,
    gelex::GenotypeMethod method,
    ChunkStats& stats) -> void
{
    const gelex::LociEncoding encoding{gelex::encode_inplace<double>(
        chunk, mode, method, 1e-12, global_start)};

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

struct GenotypeReader::EncodedChunkOutput
{
    struct Memory
    {
        Eigen::MatrixXd& data;
    };

    struct Mmap
    {
        gelex::io::BinaryWriter& writer;
        gelex::io::SectionHandle<double> genotype_handle;
    };

    explicit EncodedChunkOutput(Eigen::MatrixXd& data) : target{Memory{data}} {}

    EncodedChunkOutput(
        gelex::io::BinaryWriter& writer,
        gelex::io::SectionHandle<double> genotype_handle)
        : target{Mmap{writer, genotype_handle}}
    {
    }

    std::variant<Memory, Mmap> target;
    ChunkStats stats;
};

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

auto GenotypeReader::read_encoded_chunks(
    EncodedChunkOutput& output,
    gelex::GeneticMode mode,
    gelex::GenotypeMethod method,
    std::size_t chunk_size) -> void
{
    output.stats.means.resize(num_variants_);
    output.stats.vars.resize(num_variants_);
    output.stats.A1freqs.resize(num_variants_);
    output.stats.mono_indices.reserve(num_variants_ / 100);

    int64_t processed = 0;
    for (int64_t start = 0; start < num_variants_;)
    {
        const int64_t end
            = std::min(static_cast<int64_t>(start + chunk_size), num_variants_);
        std::visit(
            [&](auto& target)
            {
                using Target = std::decay_t<decltype(target)>;
                if constexpr (
                    std::is_same_v<Target, EncodedChunkOutput::Memory>)
                {
                    auto chunk = target.data.middleCols(start, end - start);
                    bed_.read_into<double>(chunk, start);
                    process_chunk(chunk, start, mode, method, output.stats);
                }
                else
                {
                    auto chunk = bed_.read<double>(start, end);
                    process_chunk(chunk, start, mode, method, output.stats);
                    target.writer.write(target.genotype_handle, chunk);
                }
            },
            output.target);
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

    std::ranges::sort(output.stats.mono_indices);
}

auto GenotypeReader::read_in_memory(
    gelex::GeneticMode mode,
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

    EncodedChunkOutput output{owned.data};
    read_encoded_chunks(output, mode, method, chunk_size);
    owned.mono_indices = std::move(output.stats.mono_indices);
    owned.mean = std::move(output.stats.means);
    owned.var = std::move(output.stats.vars);
    owned.A1freq = std::move(output.stats.A1freqs);

    return Genotype(std::move(owned));
}

auto GenotypeReader::read_mmap(
    gelex::GeneticMode mode,
    gelex::GenotypeMethod method,
    const std::filesystem::path& output_prefix,
    std::size_t chunk_size) -> Genotype
{
    const auto effect = gelex::EffectType::from_genetic(mode);

    auto gbin_path = output_prefix;
    gbin_path += ".gbin";

    if (std::filesystem::exists(gbin_path))
    {
        throw gelex::GelexException(
            fmt::format(
                "Output file already exists: [{}]", gbin_path.string()));
    }

    {
        gelex::io::BinaryWriter writer(gbin_path.string());
        auto genotype_handle = writer.reserve<double>(
            fmt::format("{}/genotype", effect), sample_size_, num_variants_);

        EncodedChunkOutput output{writer, genotype_handle};
        read_encoded_chunks(output, mode, method, chunk_size);

        auto stats_handle = writer.reserve<double>(
            fmt::format("{}/loci_stats", effect), num_variants_, 3);
        writer.write(stats_handle, output.stats.means);
        writer.write(stats_handle, output.stats.vars);
        writer.write(stats_handle, output.stats.A1freqs);

        if (!output.stats.mono_indices.empty())
        {
            auto mono_handle = writer.reserve<int64_t>(
                fmt::format("{}/mono_indices", effect),
                output.stats.mono_indices.size(),
                1);
            writer.write(mono_handle, output.stats.mono_indices);
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

}  // namespace gelex::genotype
