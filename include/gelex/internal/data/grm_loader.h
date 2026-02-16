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

#ifndef GELEX_INTERNAL_DATA_GRM_LOADER_H_
#define GELEX_INTERNAL_DATA_GRM_LOADER_H_

#include <cstddef>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include <mio.h>
#include <Eigen/Core>

#include "gelex/types/freq_effect.h"

namespace gelex::detail
{

class GrmLoader
{
   public:
    explicit GrmLoader(const std::filesystem::path& prefix);

    GrmLoader(const GrmLoader&) = delete;
    GrmLoader(GrmLoader&&) noexcept = default;
    GrmLoader& operator=(const GrmLoader&) = delete;
    GrmLoader& operator=(GrmLoader&&) noexcept = default;
    ~GrmLoader() = default;

    [[nodiscard]] auto load() const -> Eigen::MatrixXd;
    [[nodiscard]] auto load_unnormalized() const -> Eigen::MatrixXd;

    [[nodiscard]] auto load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> Eigen::MatrixXd;

    [[nodiscard]] auto load_unnormalized(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> Eigen::MatrixXd;

    auto load_unnormalized(
        const std::unordered_map<std::string, Eigen::Index>& id_map,
        Eigen::MatrixXd& target) const -> void;

    [[nodiscard]] auto sample_ids() const noexcept
        -> const std::vector<std::string>&
    {
        return sample_ids_;
    }

    [[nodiscard]] auto num_samples() const noexcept -> Eigen::Index
    {
        return num_samples_;
    }

    [[nodiscard]] auto type() const noexcept -> freq::GrmType { return type_; }

   private:
    std::filesystem::path bin_path_;
    std::filesystem::path id_path_;
    mio::mmap_source mmap_;
    std::vector<std::string> sample_ids_;
    Eigen::Index num_samples_{};
    freq::GrmType type_;

    auto load_sample_ids() -> void;
    auto init_mmap() -> void;

    [[nodiscard]] static auto lower_triangle_index(
        Eigen::Index i,
        Eigen::Index j) noexcept -> size_t
    {
        return (i * (i + 1) / 2) + j;
    }
};

}  // namespace gelex::detail

#endif  // GELEX_INTERNAL_DATA_GRM_LOADER_H_
