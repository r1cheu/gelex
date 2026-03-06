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

#ifndef GELEX_PIPELINE_PHENO_PIPE_H_
#define GELEX_PIPELINE_PHENO_PIPE_H_

#include <filesystem>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "gelex/data/frame/dataframe.h"
#include "gelex/data/genotype/sample_manager.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/types/fixed_effects.h"

namespace gelex
{

namespace detail
{
enum class TransformType : uint8_t
{
    None,
    DINT,
    IINT
};
}  // namespace detail

class PhenoPipe
{
   public:
    struct Config
    {
        std::filesystem::path phenotype_path;
        int phenotype_column = 3;

        std::filesystem::path bed_path;
        std::optional<std::filesystem::path> quantitative_covariates_path;
        std::optional<std::filesystem::path> discrete_covariates_path;

        detail::TransformType transform_type = detail::TransformType::None;
        double int_offset = 3.0 / 8.0;
    };

    explicit PhenoPipe(const Config& config, DataPipeObserver observer = {});
    PhenoPipe(const PhenoPipe&) = delete;
    PhenoPipe(PhenoPipe&&) noexcept = default;
    PhenoPipe& operator=(const PhenoPipe&) = delete;
    PhenoPipe& operator=(PhenoPipe&&) noexcept = default;
    ~PhenoPipe() = default;

    auto load(const std::vector<std::span<const std::string>>& extra_ids = {})
        -> void;

    auto sample_manager() -> std::shared_ptr<SampleManager>
    {
        return sample_manager_;
    }

    auto num_genotype_samples() const -> size_t
    {
        return num_genotype_samples_;
    }

    auto take_phenotype() && -> Eigen::VectorXd
    {
        return std::move(phenotype_);
    }

    auto take_fixed_effects() && -> FixedEffect
    {
        return std::move(fixed_effects_);
    }

    auto fixed_effects() const -> const FixedEffect& { return fixed_effects_; }

   private:
    auto load_phenotypes() -> void;
    auto load_covariates() -> void;
    auto intersect_samples(
        const std::vector<std::span<const std::string>>& extra_ids) -> void;
    auto finalize() -> void;
    auto apply_phenotype_transform(detail::TransformType type, double offset)
        -> void;

    Config config_;
    size_t num_genotype_samples_{};

    DataFrame<double> phenotype_frame_;
    std::optional<DataFrame<double>> qcovar_frame_;
    std::optional<DataFrame<std::string>> dcovar_frame_;
    std::string phenotype_name_;

    Eigen::VectorXd phenotype_;
    FixedEffect fixed_effects_;

    std::shared_ptr<SampleManager> sample_manager_;
    DataPipeObserver observer_;
};

}  // namespace gelex

#endif  // GELEX_PIPELINE_PHENO_PIPE_H_
