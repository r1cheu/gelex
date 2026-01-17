#ifndef GELEX_DATA_DATA_PIPE_H_
#define GELEX_DATA_DATA_PIPE_H_

#include <filesystem>
#include <memory>
#include <string>
#include <variant>
#include <vector>

#include <Eigen/Dense>

#include "../src/types/fixed_effects.h"
#include "../src/types/freq_effect.h"
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
        bool iid_only = false;
        std::string output_prefix;
        std::vector<std::filesystem::path> grm_paths;
    };

    explicit DataPipe(const Config& config);
    DataPipe(const DataPipe&) = delete;
    DataPipe(DataPipe&&) noexcept;
    DataPipe& operator=(const DataPipe&) = delete;
    DataPipe& operator=(DataPipe&&) noexcept;
    ~DataPipe();

    PhenoStats load_phenotypes();
    CovarStats load_covariates();
    std::vector<GrmStats> load_grms();
    GenotypeStats load_additive_matrix();
    GenotypeStats load_dominance_matrix();
    IntersectionStats intersect_samples();

    void finalize();

    size_t num_genotype_samples() const { return num_genotype_samples_; }

    std::shared_ptr<SampleManager> sample_manager() { return sample_manager_; }

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

    auto take_grms() && -> std::vector<gelex::freq::GeneticEffect>
    {
        return std::move(grms_);
    }

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

    std::vector<detail::GrmLoader> grm_loaders_;
    std::vector<gelex::freq::GeneticEffect> grms_;

    std::vector<std::string> fixed_effect_names_;
};

}  // namespace gelex

#endif  // GELEX_DATA_DATA_PIPE_H_
