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

#include "gelex/data/genotype_sink.h"

#include <Eigen/Core>
#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fmt/format.h>
#include <memory>
#include <new>
#include <span>
#include <utility>

#include "gelex/data/genotype.h"
#include "gelex/data/snp_stats.h"
#include "gelex/exception.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/binary_writer.h"
#include "gelex/io/snpstats.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

auto load_genotype(
    const std::filesystem::path& geno_path,
    gelex::GeneticMode mode) -> Genotype
{
    MmappedStorage mapped;
    mapped.reader = std::make_unique<gelex::BinaryReader>(geno_path.string());

    auto geno_map
        = mapped.reader->to_map<double>(fmt::format("{}/genotype", mode));
    new (&mapped.view) MmappedStorage::MapType(
        geno_map.data(), geno_map.rows(), geno_map.cols());

    mapped.stats = read_snp_stats(*mapped.reader, mode);
    return Genotype(std::move(mapped));
}

InMemorySink::InMemorySink(
    std::span<const gelex::GeneticMode> modes,
    Eigen::Index sample_size,
    Eigen::Index num_variants)
    : modes_(modes.begin(), modes.end())
{
    const std::size_t nm = modes_.size();
    owned_.resize(nm);
    targets_.reserve(nm);
    try
    {
        for (auto& storage : owned_)
        {
            storage.data.resize(sample_size, num_variants);
        }
    }
    catch (const std::bad_alloc&)
    {
        const double gb = static_cast<double>(nm)
                          * static_cast<double>(sample_size)
                          * static_cast<double>(num_variants) * sizeof(double)
                          / 1024.0 / 1024.0 / 1024.0;
        throw gelex::GelexException(
            fmt::format(
                "Memory allocation failed for {} genotype matrix/matrices "
                "({} x {}). Requires approx {:.2f} GB RAM.",
                nm,
                sample_size,
                num_variants,
                gb));
    }
}

auto InMemorySink::chunk_targets(Eigen::Index start, Eigen::Index cols)
    -> std::span<Eigen::Ref<Eigen::MatrixXd>>
{
    targets_.clear();
    for (auto& storage : owned_)
    {
        targets_.emplace_back(storage.data.middleCols(start, cols));
    }
    return std::span{targets_};
}

auto InMemorySink::commit_chunk() -> void {}

auto InMemorySink::finalize(std::span<gelex::SnpStats> stats)
    -> gelex::ModeMap<gelex::Genotype>
{
    gelex::ModeMap<gelex::Genotype> result;
    for (std::size_t k = 0; k < modes_.size(); ++k)
    {
        owned_[k].stats = std::move(stats[k]);
        result.emplace(modes_[k], Genotype(std::move(owned_[k])));
    }
    return result;
}

MmapSink::MmapSink(
    std::span<const gelex::GeneticMode> modes,
    Eigen::Index sample_size,
    Eigen::Index num_variants,
    const std::filesystem::path& output_prefix)
    : modes_(modes.begin(), modes.end()), sample_size_(sample_size)
{
    const std::size_t nm = modes_.size();
    paths_.resize(nm);
    writers_.reserve(nm);
    handles_.reserve(nm);
    for (std::size_t k = 0; k < nm; ++k)
    {
        auto geno_path = output_prefix;
        geno_path += fmt::format(".{}.geno", modes_[k]);
        if (std::filesystem::exists(geno_path))
        {
            throw gelex::GelexException(
                fmt::format(
                    "Output file already exists: [{}]", geno_path.string()));
        }
        paths_[k] = geno_path;
        writers_.push_back(
            std::make_unique<gelex::BinaryWriter>(geno_path.string()));
        handles_.push_back(
            writers_[k]->reserve<double>(
                fmt::format("{}/genotype", modes_[k]),
                sample_size,
                num_variants));
    }
}

auto MmapSink::chunk_targets(Eigen::Index /*start*/, Eigen::Index cols)
    -> std::span<Eigen::Ref<Eigen::MatrixXd>>
{
    temps_.clear();
    temps_.reserve(modes_.size());
    for (std::size_t k = 0; k < modes_.size(); ++k)
    {
        temps_.emplace_back(sample_size_, cols);
    }
    targets_.clear();
    targets_.reserve(temps_.size());
    for (auto& temp : temps_)
    {
        targets_.emplace_back(temp);
    }
    return std::span{targets_};
}

auto MmapSink::commit_chunk() -> void
{
    for (std::size_t k = 0; k < writers_.size(); ++k)
    {
        writers_[k]->write(handles_[k], temps_[k]);
    }
}

auto MmapSink::finalize(std::span<gelex::SnpStats> stats)
    -> gelex::ModeMap<gelex::Genotype>
{
    const std::size_t nm = modes_.size();
    for (std::size_t k = 0; k < nm; ++k)
    {
        write_snp_stats(*writers_[k], modes_[k], stats[k]);
    }
    writers_.clear();  // flush and close every file before mapping it back

    gelex::ModeMap<gelex::Genotype> result;
    for (std::size_t k = 0; k < nm; ++k)
    {
        result.emplace(modes_[k], load_genotype(paths_[k], modes_[k]));
    }
    return result;
}

}  // namespace gelex
