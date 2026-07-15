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

#ifndef GELEX_DATA_GENOTYPE_SINK_H_
#define GELEX_DATA_GENOTYPE_SINK_H_

#include <Eigen/Core>
#include <filesystem>
#include <memory>
#include <span>
#include <vector>

#include "gelex/data/genotype.h"
#include "gelex/data/snp_stats.h"
#include "gelex/io/binary_writer.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

// Where an encoding pass lands each chunk's columns and how it assembles the
// finished Genotype set. This is the axis that separates resident matrices from
// on-disk mmap; the pass driving it is shared and mode-agnostic.
class GenotypeSink
{
   public:
    GenotypeSink() = default;
    GenotypeSink(const GenotypeSink&) = delete;
    GenotypeSink(GenotypeSink&&) = delete;
    auto operator=(const GenotypeSink&) -> GenotypeSink& = delete;
    auto operator=(GenotypeSink&&) -> GenotypeSink& = delete;
    virtual ~GenotypeSink() = default;

    // Per-mode write targets for the columns [start, start + cols). The
    // returned Refs stay valid until the next chunk_targets or commit_chunk
    // call.
    virtual auto chunk_targets(Eigen::Index start, Eigen::Index cols)
        -> std::span<Eigen::Ref<Eigen::MatrixXd>> = 0;

    // Persists the just-encoded chunk; a no-op for sinks that encode in place.
    virtual auto commit_chunk() -> void = 0;

    virtual auto finalize(std::span<gelex::SnpStats> stats)
        -> gelex::ModeMap<gelex::Genotype> = 0;
};

class InMemorySink final : public GenotypeSink
{
   public:
    InMemorySink(
        std::span<const gelex::GeneticMode> modes,
        Eigen::Index sample_size,
        Eigen::Index num_variants);

    auto chunk_targets(Eigen::Index start, Eigen::Index cols)
        -> std::span<Eigen::Ref<Eigen::MatrixXd>> override;
    auto commit_chunk() -> void override;
    auto finalize(std::span<gelex::SnpStats> stats)
        -> gelex::ModeMap<gelex::Genotype> override;

   private:
    std::vector<gelex::GeneticMode> modes_;
    std::vector<gelex::OwnedStorage> owned_;
    std::vector<Eigen::Ref<Eigen::MatrixXd>> targets_;
};

class MmapSink final : public GenotypeSink
{
   public:
    MmapSink(
        std::span<const gelex::GeneticMode> modes,
        Eigen::Index sample_size,
        Eigen::Index num_variants,
        const std::filesystem::path& output_prefix);

    auto chunk_targets(Eigen::Index start, Eigen::Index cols)
        -> std::span<Eigen::Ref<Eigen::MatrixXd>> override;
    auto commit_chunk() -> void override;
    auto finalize(std::span<gelex::SnpStats> stats)
        -> gelex::ModeMap<gelex::Genotype> override;

   private:
    std::vector<gelex::GeneticMode> modes_;
    Eigen::Index sample_size_;
    std::vector<std::filesystem::path> paths_;
    std::vector<std::unique_ptr<gelex::BinaryWriter>> writers_;
    std::vector<gelex::SectionHandle<double>> handles_;
    std::vector<Eigen::MatrixXd> temps_;
    std::vector<Eigen::Ref<Eigen::MatrixXd>> targets_;
};

}  // namespace gelex

#endif  // GELEX_DATA_GENOTYPE_SINK_H_
