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
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "gelex/data/encode/encoder.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"
#include "gelex/data/genotype.h"
#include "gelex/data/genotype_method.h"
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

// Fuses each variant straight from its packed form into the chunk column,
// filling the shared SnpStats in one pass instead of decoding to dosage first.
auto encode_chunk(
    const gelex::LocusEncoder& encoder,
    Eigen::Ref<Eigen::MatrixXd> chunk,
    Eigen::Index global_start,
    const gelex::EncodingSpec& spec,
    gelex::SnpStats& stats) -> void
{
    for (Eigen::Index c = 0; c < chunk.cols(); ++c)
    {
        const Eigen::Index marker = global_start + c;
        const gelex::LocusStats locus_stats{encoder.count(marker)};
        const gelex::LocusEncoding encoding{
            encoder.encoding(marker, locus_stats, spec)};
        encoder.expand(marker, encoding, chunk.col(c));

        stats.code.col(marker) = encoding.code.matrix();
        stats.A1freq[marker]
            = locus_stats.has_nonmissing() ? locus_stats.A1freq() : 0.0;

        if (encoding.valid)
        {
            stats.valid_indices.push_back(static_cast<int64_t>(marker));
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

struct GenotypeReader::EncodedChunkOutput
{
    struct Memory
    {
        Eigen::MatrixXd& data;
    };

    struct Mmap
    {
        gelex::BinaryWriter& writer;
        gelex::SectionHandle<double> genotype_handle;
    };

    explicit EncodedChunkOutput(Eigen::MatrixXd& data) : target{Memory{data}} {}

    EncodedChunkOutput(
        gelex::BinaryWriter& writer,
        gelex::SectionHandle<double> genotype_handle)
        : target{Mmap{writer, genotype_handle}}
    {
    }

    std::variant<Memory, Mmap> target;
    gelex::SnpStats stats;
};

GenotypeReader::GenotypeReader(
    const gelex::Bed& bed,
    gelex::GenoObserver observer)
    : bed_(bed),
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
    output.stats.code.resize(3, num_variants_);
    output.stats.A1freq.resize(num_variants_);
    output.stats.valid_indices.reserve(num_variants_);

    const gelex::LocusEncoder encoder{bed_};
    const gelex::EncodingSpec spec{encoding_spec_from_method(mode, method)};

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
                    encode_chunk(encoder, chunk, start, spec, output.stats);
                }
                else
                {
                    Eigen::MatrixXd chunk(sample_size_, end - start);
                    encode_chunk(encoder, chunk, start, spec, output.stats);
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

    std::ranges::sort(output.stats.valid_indices);
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
    owned.stats = std::move(output.stats);

    return Genotype(std::move(owned));
}

auto GenotypeReader::read_mmap(
    gelex::GeneticMode mode,
    gelex::GenotypeMethod method,
    const std::filesystem::path& output_prefix,
    std::size_t chunk_size) -> Genotype
{
    auto geno_path = output_prefix;
    geno_path += ".geno";

    if (std::filesystem::exists(geno_path))
    {
        throw gelex::GelexException(
            fmt::format(
                "Output file already exists: [{}]", geno_path.string()));
    }

    {
        gelex::BinaryWriter writer(geno_path.string());
        auto genotype_handle = writer.reserve<double>(
            fmt::format("{}/genotype", mode), sample_size_, num_variants_);

        EncodedChunkOutput output{writer, genotype_handle};
        read_encoded_chunks(output, mode, method, chunk_size);

        write_snp_stats(writer, mode, output.stats);
    }

    return Genotype(load_mmapped(geno_path, mode));
}

auto GenotypeReader::read(
    const std::filesystem::path& geno_path,
    gelex::GeneticMode mode) -> Genotype
{
    return Genotype(load_mmapped(geno_path, mode));
}

}  // namespace gelex
