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

#ifndef GELEX_DATA_GRM_H_
#define GELEX_DATA_GRM_H_

#include <span>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/bed.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/locus_encoding.h"
#include "gelex/infra/logging/grm_event.h"
#include "gelex/infra/logging/notify.h"

namespace gelex
{

struct GrmResult
{
    Eigen::MatrixXd grm;
    double denominator;
};

class GRM
{
   public:
    explicit GRM(const std::string& bfile_prefix);
    GRM(const GRM&) = delete;
    GRM(GRM&&) noexcept = default;
    GRM& operator=(const GRM&) = delete;
    GRM& operator=(GRM&&) noexcept = default;

    ~GRM() = default;

    template <GeneticMode GT>
    auto compute(
        GenotypeMethod method,
        Eigen::Index chunk_size,
        const GrmObserver& observer = {}) -> GrmResult;

    template <GeneticMode GT>
    auto compute(
        GenotypeMethod method,
        const std::vector<std::pair<Eigen::Index, Eigen::Index>>& ranges,
        Eigen::Index chunk_size,
        const GrmObserver& observer = {}) -> GrmResult;

    [[nodiscard]] auto sample_ids() const -> std::span<const std::string>
    {
        return sample_index_.keys();
    }

    [[nodiscard]] auto num_snps() const -> Eigen::Index
    {
        return bed_.num_snps();
    }

   private:
    dataframe::Index<std::string> sample_index_;
    Bed bed_;

    static auto update_grm(
        Eigen::Ref<Eigen::MatrixXd> grm,
        const Eigen::Ref<const Eigen::MatrixXd>& genotype) -> void;
};

template <GeneticMode GT>
auto GRM::compute(
    GenotypeMethod method,
    Eigen::Index chunk_size,
    const GrmObserver& observer) -> GrmResult
{
    return compute<GT>(method, {{0, bed_.num_snps()}}, chunk_size, observer);
}

template <GeneticMode GT>
auto GRM::compute(
    GenotypeMethod method,
    const std::vector<std::pair<Eigen::Index, Eigen::Index>>& ranges,
    Eigen::Index chunk_size,
    const GrmObserver& observer) -> GrmResult
{
    const Eigen::Index n = bed_.num_samples();
    Eigen::MatrixXd grm = Eigen::MatrixXd::Zero(n, n);

    Eigen::Index total_snps_to_process = 0;
    for (const auto& [start, end] : ranges)
    {
        total_snps_to_process += (end - start);
    }

    Eigen::Index processed_snps = 0;
    for (const auto& [range_start, range_end] : ranges)
    {
        for (Eigen::Index start_col = range_start; start_col < range_end;
             start_col += chunk_size)
        {
            const Eigen::Index end_col
                = std::min(start_col + chunk_size, range_end);
            Eigen::MatrixXd genotype_chunk
                = bed_.read<double>(start_col, end_col);

            encode_inplace<double>(genotype_chunk, GT, method);
            update_grm(grm, genotype_chunk);

            processed_snps += (end_col - start_col);
            notify(
                observer,
                GrmProgressEvent{
                    static_cast<size_t>(processed_snps),
                    static_cast<size_t>(total_snps_to_process),
                    false});
        }
    }

    double denominator = grm.trace() / static_cast<double>(n);

    return {grm, denominator};
}
}  // namespace gelex

#endif  // GELEX_DATA_GRM_H_
