#pragma once

#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <variant>
#include <vector>

#include <Eigen/Dense>

#include "../src/types/fixed_effects.h"
#include "gelex/data/genotype_loader.h"
#include "gelex/data/genotype_matrix.h"
#include "gelex/data/genotype_mmap.h"
#include "gelex/data/genotype_pipe.h"
#include "gelex/data/sample_manager.h"

namespace gelex
{
namespace detail
{
class PhenotypeLoader;
class QuantitativeCovariateLoader;
class DiscreteCovariateLoader;
class GrmLoader;
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
        bool iid_only = false;
        std::string output_prefix;
        std::filesystem::path additive_grm_path;   // GRM file prefix (加性)
        std::filesystem::path dominance_grm_path;  // GRM file prefix (显性)
    };

    explicit DataPipe(const Config& config);
    DataPipe(const DataPipe&) = delete;
    DataPipe(DataPipe&&) noexcept;
    DataPipe& operator=(const DataPipe&) = delete;
    DataPipe& operator=(DataPipe&&) noexcept;
    ~DataPipe();

    PhenoStats load_phenotypes();
    CovarStats load_covariates();
    GrmStats load_additive_grm();
    GrmStats load_dominance_grm();
    GenotypeStats load_additive_matrix();
    GenotypeStats load_dominance_matrix();
    IntersectionStats intersect_samples();

    void finalize();

    size_t num_genotype_samples() const { return num_genotype_samples_; }

    Eigen::VectorXd take_phenotype() && { return std::move(phenotype_); }
    FixedEffect take_fixed_effects() && { return std::move(fixed_effects_); }
    const FixedEffect& fixed_effects() const { return fixed_effects_; }
    std::variant<GenotypeMap, GenotypeMatrix> take_additive_matrix() &&
    {
        return std::move(*additive_matrix_);
    }
    std::variant<GenotypeMap, GenotypeMatrix> take_dominance_matrix() &&
    {
        return std::move(*dominance_matrix_);
    }
    bool has_dominance_matrix() const { return dominance_matrix_ != nullptr; }

    Eigen::MatrixXd take_additive_grm() && { return std::move(*additive_grm_); }
    Eigen::MatrixXd take_dominance_grm() &&
    {
        return std::move(*dominance_grm_);
    }
    bool has_additive_grm() const { return additive_grm_ != nullptr; }
    bool has_dominance_grm() const { return dominance_grm_ != nullptr; }

    const std::vector<std::string>& fixed_effect_names() const;

   private:
    DataPipe() = default;

    template <VariantProcessor Processor, typename TargetPtr>
    void load_genotype_impl(const std::string& suffix, TargetPtr& target)
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

    std::unique_ptr<detail::PhenotypeLoader> phenotype_loader_;
    std::unique_ptr<detail::QuantitativeCovariateLoader> qcovar_loader_;
    std::unique_ptr<detail::DiscreteCovariateLoader> dcovar_loader_;

    Eigen::VectorXd phenotype_;
    FixedEffect fixed_effects_;

    std::shared_ptr<SampleManager> sample_manager_;

    std::unique_ptr<std::variant<GenotypeMap, GenotypeMatrix>> additive_matrix_;
    std::unique_ptr<std::variant<GenotypeMap, GenotypeMatrix>>
        dominance_matrix_;

    std::unique_ptr<detail::GrmLoader> additive_grm_loader_;
    std::unique_ptr<detail::GrmLoader> dominance_grm_loader_;
    std::unique_ptr<Eigen::MatrixXd> additive_grm_;
    std::unique_ptr<Eigen::MatrixXd> dominance_grm_;

    std::vector<std::string> fixed_effect_names_;
};

}  // namespace gelex
