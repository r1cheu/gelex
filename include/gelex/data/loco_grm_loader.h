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

#ifndef GELEX_DATA_LOCO_GRM_LOADER_H_
#define GELEX_DATA_LOCO_GRM_LOADER_H_

#include <filesystem>
#include <string>
#include <unordered_map>

#include <Eigen/Core>

namespace gelex
{
namespace detail
{
class GrmLoader;
}

class LocoGRMLoader
{
   public:
    LocoGRMLoader(
        const std::filesystem::path& whole_grm_prefix,
        const std::unordered_map<std::string, Eigen::Index>& id_map);

    LocoGRMLoader(const LocoGRMLoader&) = delete;
    LocoGRMLoader(LocoGRMLoader&&) noexcept = default;
    LocoGRMLoader& operator=(const LocoGRMLoader&) = delete;
    LocoGRMLoader& operator=(LocoGRMLoader&&) noexcept = default;
    ~LocoGRMLoader();

    /**
     * @brief Load the LOCO GRM for a specific chromosome.
     *
     * Formula: G_loco = (G_whole - G_i) / (K_whole - K_i)
     *
     * @param chr_grm_prefix Path prefix for the chromosome-specific GRM files.
     * @param id_map Map for sample IDs to matrix indices.
     * @param target Output matrix to store the calculated LOCO GRM.
     */
    auto load_loco_grm(
        const std::filesystem::path& chr_grm_prefix,
        const std::unordered_map<std::string, Eigen::Index>& id_map,
        Eigen::MatrixXd& target) const -> void;

    [[nodiscard]] auto load_loco_grm(
        const std::filesystem::path& chr_grm_prefix,
        const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> Eigen::MatrixXd;

    [[nodiscard]] auto num_samples() const noexcept -> Eigen::Index;

   private:
    Eigen::MatrixXd
        g_whole_;           // Pre-filtered unnormalized whole matrix (Z * Z')
    double k_whole_{};      // K_whole (= trace(g_whole_) / n)
    double trace_whole_{};  // trace(g_whole_), saved for LOCO calculation
    mutable Eigen::MatrixXd g_chr_buffer_;  // Buffer for chromosome GRM
};

}  // namespace gelex

#endif  // GELEX_DATA_LOCO_GRM_LOADER_H_
