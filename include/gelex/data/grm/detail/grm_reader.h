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

#ifndef GELEX_DATA_GRM_DETAIL_GRM_READER_H_
#define GELEX_DATA_GRM_DETAIL_GRM_READER_H_

#include <cstddef>
#include <filesystem>
#include <string>

#include <mio.h>
#include <Eigen/Core>

#include "gelex/data/dataframe/index.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::grm::detail
{

class GrmReader
{
   public:
    explicit GrmReader(const std::filesystem::path& prefix);

    GrmReader(const GrmReader&) = delete;
    GrmReader(GrmReader&&) noexcept = default;
    GrmReader& operator=(const GrmReader&) = delete;
    GrmReader& operator=(GrmReader&&) noexcept = default;
    ~GrmReader() = default;

    [[nodiscard]] auto load() const -> Eigen::MatrixXd;
    [[nodiscard]] auto load_unnormalized() const -> Eigen::MatrixXd;

    [[nodiscard]] auto load(const dataframe::Index<std::string>& sample_index)
        const -> Eigen::MatrixXd;

    [[nodiscard]] auto load_unnormalized(
        const dataframe::Index<std::string>& sample_index) const
        -> Eigen::MatrixXd;

    auto load_unnormalized(
        const dataframe::Index<std::string>& sample_index,
        Eigen::MatrixXd& target) const -> void;

    [[nodiscard]] auto sample_index() const noexcept
        -> const dataframe::Index<std::string>&
    {
        return sample_index_;
    }

    [[nodiscard]] auto num_samples() const noexcept -> Eigen::Index
    {
        return static_cast<Eigen::Index>(sample_index_.size());
    }

    [[nodiscard]] auto type() const noexcept -> GeneticMode { return type_; }

   private:
    std::filesystem::path bin_path_;
    mio::mmap_source mmap_;
    dataframe::Index<std::string> sample_index_;
    GeneticMode type_;

    auto init_mmap() -> void;

    [[nodiscard]] static auto lower_triangle_index(
        Eigen::Index i,
        Eigen::Index j) noexcept -> size_t
    {
        return (i * (i + 1) / 2) + j;
    }
};

}  // namespace gelex::grm::detail

#endif  // GELEX_DATA_GRM_DETAIL_GRM_READER_H_
