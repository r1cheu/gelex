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

#ifndef GELEX_DATA_GENOTYPE_H_
#define GELEX_DATA_GENOTYPE_H_

#include <Eigen/Core>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <utility>
#include <variant>
#include <vector>

#include "gelex/data/snp_stats.h"
#include "gelex/io/binary_reader.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::test
{
class GenotypeBuilder;
}  // namespace gelex::test

namespace gelex
{

class Genotype;
class InMemorySink;

auto load_genotype(
    const std::filesystem::path& geno_path,
    gelex::GeneticMode mode) -> Genotype;

struct OwnedStorage
{
    Eigen::MatrixXd data;
    SnpStats stats;

    OwnedStorage() = default;
    OwnedStorage(const OwnedStorage&) = delete;
    OwnedStorage(OwnedStorage&&) noexcept = default;
    auto operator=(const OwnedStorage&) -> OwnedStorage& = delete;
    auto operator=(OwnedStorage&&) noexcept -> OwnedStorage& = default;
    ~OwnedStorage() = default;
};

struct MmappedStorage
{
    using MapType = Eigen::Map<const Eigen::MatrixXd, Eigen::Unaligned>;

    std::unique_ptr<gelex::BinaryReader> reader;
    MapType view{nullptr, 0, 0};
    SnpStats stats;

    MmappedStorage() = default;
    MmappedStorage(const MmappedStorage&) = delete;
    MmappedStorage(MmappedStorage&&) noexcept = default;
    auto operator=(const MmappedStorage&) -> MmappedStorage& = delete;
    auto operator=(MmappedStorage&&) noexcept -> MmappedStorage& = default;
    ~MmappedStorage() = default;
};

}  // namespace gelex

namespace gelex
{

class Genotype
{
   public:
    Genotype(const Genotype&) = delete;
    Genotype(Genotype&&) noexcept = default;
    auto operator=(const Genotype&) -> Genotype& = delete;
    auto operator=(Genotype&&) noexcept -> Genotype& = default;
    ~Genotype() = default;

    [[nodiscard]] auto matrix() const noexcept
        -> Eigen::Ref<const Eigen::MatrixXd>;

    [[nodiscard]] auto stats() const noexcept -> const SnpStats&;

    [[nodiscard]] auto A1freq() const noexcept -> const Eigen::VectorXd&;

    [[nodiscard]] auto valid_indices() const noexcept
        -> const std::vector<int64_t>&;

    [[nodiscard]] auto num_invalid() const noexcept -> int64_t;

    [[nodiscard]] auto rows() const noexcept -> int64_t;

    [[nodiscard]] auto cols() const noexcept -> int64_t;

   private:
    friend class InMemorySink;
    friend auto load_genotype(
        const std::filesystem::path& geno_path,
        gelex::GeneticMode mode) -> Genotype;
    friend class ::gelex::test::GenotypeBuilder;

    explicit Genotype(OwnedStorage&& owned) noexcept
        : storage_(std::move(owned))
    {
    }

    explicit Genotype(MmappedStorage&& mapped) noexcept
        : storage_(std::move(mapped))
    {
    }

    std::variant<OwnedStorage, MmappedStorage> storage_;
};

}  // namespace gelex

#endif  // GELEX_DATA_GENOTYPE_H_
