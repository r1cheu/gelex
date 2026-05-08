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

#ifndef GELEX_DATA_PIPE_PHENO_H_
#define GELEX_DATA_PIPE_PHENO_H_

#include <filesystem>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/infra/logging/pheno_event.h"
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

        std::optional<std::filesystem::path> quantitative_covariates_path;
        std::optional<std::filesystem::path> discrete_covariates_path;

        detail::TransformType transform_type = detail::TransformType::None;
        double int_offset = 3.0 / 8.0;
    };

    explicit PhenoPipe(const Config& config, PhenoObserver observer = {});
    PhenoPipe(const PhenoPipe&) = delete;
    PhenoPipe(PhenoPipe&&) noexcept = default;
    PhenoPipe& operator=(const PhenoPipe&) = delete;
    PhenoPipe& operator=(PhenoPipe&&) noexcept = default;
    ~PhenoPipe() = default;

    auto sample_indices() const
        -> std::vector<const dataframe::Index<std::string>*>;

    auto load(const dataframe::Index<std::string>& common) -> void;

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
    auto apply_phenotype_transform(detail::TransformType type, double offset)
        -> void;

    Config config_;

    std::optional<dataframe::DataFrame<std::string>> phenotype_frame_;
    std::optional<dataframe::DataFrame<std::string>> qcovar_frame_;
    std::optional<dataframe::DataFrame<std::string>> dcovar_frame_;

    Eigen::VectorXd phenotype_;
    FixedEffect fixed_effects_;

    PhenoObserver observer_;
};

}  // namespace gelex

#endif  // GELEX_DATA_PIPE_PHENO_H_
