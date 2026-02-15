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

#ifndef GELEX_DATA_DATA_PIPE_H_
#define GELEX_DATA_DATA_PIPE_H_

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Dense>

#include "../src/types/fixed_effects.h"
#include "../src/types/freq_effect.h"
#include "gelex/data/dataframe.h"
#include "gelex/data/genotype_loader.h"
#include "gelex/data/genotype_matrix.h"
#include "gelex/data/genotype_method_dispatch.h"
#include "gelex/data/genotype_mmap.h"
#include "gelex/data/genotype_pipe.h"
#include "gelex/data/sample_manager.h"

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

struct PhenoStats
{
    size_t samples_loaded;
    std::string trait_name;
};

struct CovarStats
{
    size_t qcovar_loaded;
    size_t dcovar_loaded;
    std::vector<std::string> q_names;
    std::vector<std::string> d_names;
};

struct IntersectionStats
{
    size_t total_samples;
    size_t common_samples;
    size_t excluded_samples;
};

struct GenotypeStats
{
    int64_t num_snps;
    int64_t monomorphic_snps;
};

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
            = GenotypeProcessMethod::OrthStandardize;
    };

    explicit DataPipe(const Config& config);
    DataPipe(const DataPipe&) = delete;
    DataPipe(DataPipe&&) noexcept;
    DataPipe& operator=(const DataPipe&) = delete;
    DataPipe& operator=(DataPipe&&) noexcept;
    ~DataPipe();

    auto load_phenotypes() -> PhenoStats;
    auto load_covariates() -> CovarStats;
    auto load_grms() -> std::vector<GrmStats>;
    auto load_additive_matrix() -> GenotypeStats;
    auto load_dominance_matrix() -> GenotypeStats;
    auto intersect_samples() -> IntersectionStats;

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

    template <GenotypeProcessor Processor, typename TargetPtr>
    auto load_genotype_impl(const std::string& suffix, TargetPtr& target)
        -> void
    {
        auto assign_to_target = [&](auto&& result_value)
        {
            using VariantType = typename TargetPtr::element_type;
            target = std::make_unique<VariantType>(
                std::forward<decltype(result_value)>(result_value));
        };

        if (config_.use_mmap)
        {
            std::string file_path = config_.output_prefix + suffix;
            auto pipe = gelex::GenotypePipe(
                config_.bed_path, sample_manager_, file_path);
            auto result = pipe.template process<Processor>(config_.chunk_size);
            assign_to_target(std::move(result));
        }
        else
        {
            auto loader
                = gelex::GenotypeLoader(config_.bed_path, sample_manager_);
            auto result
                = loader.template process<Processor>(config_.chunk_size);
            assign_to_target(std::move(result));
        }
    }

    Config config_;
    size_t num_genotype_samples_;

    std::optional<DataFrame<double>> phenotype_frame_;
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

    auto apply_phenotype_transform(detail::TransformType type, double offset)
        -> void;
};

}  // namespace gelex

#endif  // GELEX_DATA_DATA_PIPE_H_
