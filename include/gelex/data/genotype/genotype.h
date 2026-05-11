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

#ifndef GELEX_DATA_GENOTYPE_GENOTYPE_H_
#define GELEX_DATA_GENOTYPE_GENOTYPE_H_

#include <cstdint>
#include <memory>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/io/detail/binary_reader.h"

namespace gelex::test
{
class GenotypeBuilder;
}  // namespace gelex::test

namespace gelex::genotype
{

class GenotypeReader;

struct OwnedStorage
{
    Eigen::MatrixXd data;
    Eigen::VectorXd mean;
    Eigen::VectorXd stddev;
    std::vector<int64_t> mono_indices;

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

    std::unique_ptr<gelex::io::detail::BinaryReader> reader;
    MapType view{nullptr, 0, 0};
    Eigen::VectorXd mean;
    Eigen::VectorXd stddev;
    std::vector<int64_t> mono_indices;

    MmappedStorage() = default;
    MmappedStorage(const MmappedStorage&) = delete;
    MmappedStorage(MmappedStorage&&) noexcept = default;
    auto operator=(const MmappedStorage&) -> MmappedStorage& = delete;
    auto operator=(MmappedStorage&&) noexcept -> MmappedStorage& = default;
    ~MmappedStorage() = default;
};

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

    [[nodiscard]] auto mean() const noexcept -> const Eigen::VectorXd&;

    [[nodiscard]] auto stddev() const noexcept -> const Eigen::VectorXd&;

    [[nodiscard]] auto mono_indices() const noexcept
        -> const std::vector<int64_t>&;

    [[nodiscard]] auto num_mono() const noexcept -> int64_t;

    [[nodiscard]] auto rows() const noexcept -> int64_t;

    [[nodiscard]] auto cols() const noexcept -> int64_t;

    [[nodiscard]] auto is_monomorphic(Eigen::Index marker_idx) const noexcept
        -> bool;

   private:
    friend class GenotypeReader;
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

}  // namespace gelex::genotype

#endif  // GELEX_DATA_GENOTYPE_GENOTYPE_H_
