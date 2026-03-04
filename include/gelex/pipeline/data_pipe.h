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

#ifndef GELEX_PIPELINE_DATA_PIPE_H_
#define GELEX_PIPELINE_DATA_PIPE_H_

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Dense>

#include "gelex/data/frame/dataframe.h"
#include "gelex/data/genotype/genotype_loader.h"
#include "gelex/data/genotype/genotype_matrix.h"
#include "gelex/data/genotype/genotype_mmap.h"
#include "gelex/data/genotype/genotype_pipe.h"
#include "gelex/data/genotype/genotype_processor.h"
#include "gelex/data/genotype/sample_manager.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/types/fixed_effects.h"
#include "gelex/types/freq_effect.h"

namespace gelex
{
namespace detail
{
class GrmLoader;

enum class TransformType : uint8_t
{
    None,
    DINT,
    IINT
};
}  // namespace detail

struct GrmStats
{
    size_t samples_in_file;
    freq::GrmType type;
};

class DataPipe
{
   public:
    struct Config
    {
        std::filesystem::path phenotype_path;
        int phenotype_column = 3;
        std::filesystem::path bed_path;
        bool use_dominance_effect = false;
        bool use_mmap = false;
        int chunk_size = 10000;
        std::filesystem::path qcovar_path;
        std::filesystem::path dcovar_path;
        std::string output_prefix;
        std::vector<std::filesystem::path> grm_paths;
        detail::TransformType transform_type = detail::TransformType::None;
        double int_offset = 3.0 / 8.0;
        GenotypeProcessMethod genotype_method
            = GenotypeProcessMethod::OrthStandardizeHWE;
    };

    explicit DataPipe(const Config& config, DataPipeObserver observer = {});
    DataPipe(const DataPipe&) = delete;
    DataPipe(DataPipe&&) noexcept;
    DataPipe& operator=(const DataPipe&) = delete;
    DataPipe& operator=(DataPipe&&) noexcept;
    ~DataPipe();

    auto load_phenotypes() -> void;
    auto load_covariates() -> void;
    auto load_grms() -> std::vector<GrmStats>;
    auto load_additive_matrix() -> void;
    auto load_dominance_matrix() -> void;
    auto intersect_samples() -> void;

    auto finalize() -> void;

    auto num_genotype_samples() const -> size_t
    {
        return num_genotype_samples_;
    }

    auto sample_manager() -> std::shared_ptr<SampleManager>
    {
        return sample_manager_;
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
    auto take_additive_matrix() && -> std::variant<GenotypeMap, GenotypeMatrix>
    {
        return std::move(*additive_matrix_);
    }
    auto take_dominance_matrix() && -> std::variant<GenotypeMap, GenotypeMatrix>
    {
        return std::move(*dominance_matrix_);
    }
    auto has_dominance_matrix() const -> bool
    {
        return dominance_matrix_ != nullptr;
    }

    auto take_grms() && -> std::vector<gelex::freq::GeneticEffect>
    {
        return std::move(grms_);
    }

   private:
    DataPipe() = default;

    using GenotypeMatrixPtr
        = std::unique_ptr<std::variant<GenotypeMap, GenotypeMatrix>>;

    template <GeneticEffectType GT>
    auto load_genotype_impl(
        const std::string& suffix,
        GenotypeProcessMethod method,
        GenotypeMatrixPtr& target) -> void
    {
        if (config_.use_mmap)
        {
            std::string file_path = config_.output_prefix + suffix;
            auto pipe = gelex::GenotypePipe(
                config_.bed_path, sample_manager_, file_path);
            target
                = std::make_unique<std::variant<GenotypeMap, GenotypeMatrix>>(
                    pipe.process<GT>(method, config_.chunk_size));
        }
        else
        {
            auto loader
                = gelex::GenotypeLoader(config_.bed_path, sample_manager_);
            target
                = std::make_unique<std::variant<GenotypeMap, GenotypeMatrix>>(
                    loader.process<GT>(method, config_.chunk_size));
        }
    }

    Config config_;
    size_t num_genotype_samples_;

    DataFrame<double> phenotype_frame_;
    std::optional<DataFrame<double>> qcovar_frame_;
    std::optional<DataFrame<std::string>> dcovar_frame_;
    std::string phenotype_name_;

    Eigen::VectorXd phenotype_;
    FixedEffect fixed_effects_;

    std::shared_ptr<SampleManager> sample_manager_;

    std::unique_ptr<std::variant<GenotypeMap, GenotypeMatrix>> additive_matrix_;
    std::unique_ptr<std::variant<GenotypeMap, GenotypeMatrix>>
        dominance_matrix_;

    std::vector<detail::GrmLoader> grm_loaders_;
    std::vector<gelex::freq::GeneticEffect> grms_;

    DataPipeObserver observer_;

    auto apply_phenotype_transform(detail::TransformType type, double offset)
        -> void;
};

}  // namespace gelex

#endif  // GELEX_PIPELINE_DATA_PIPE_H_
