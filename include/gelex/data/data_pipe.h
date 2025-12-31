#pragma once

#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <variant>
#include <vector>

#include <Eigen/Dense>

#include "Eigen/Core"

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
class QcovarLoader;
class DcovarLoader;
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
    size_t snps_loaded;
    size_t snps_total;
    bool dominance_loaded;
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
    };

    explicit DataPipe(const Config& config);
    DataPipe(const DataPipe&) = delete;
    DataPipe(DataPipe&&) noexcept;
    DataPipe& operator=(const DataPipe&) = delete;
    DataPipe& operator=(DataPipe&&) noexcept;
    ~DataPipe();

    PhenoStats load_phenotypes();
    CovarStats load_covariates();
    IntersectionStats intersect_samples();

    // build_matrices accepts a callback for progress updates.
    // The callback receives (processed_snps, total_snps) for the current matrix
    // being built.
    GenotypeStats build_matrices(
        std::function<void(size_t, size_t)> progress_callback = nullptr);

    Eigen::VectorXd take_phenotype() && { return std::move(phenotype_); }
    Eigen::MatrixXd take_fixed_effects() &&
    {
        return std::move(fixed_effects_);
    }
    std::variant<GenotypeMap, GenotypeMatrix> take_additive_matrix() &&
    {
        return std::move(*additive_matrix_);
    }
    std::variant<GenotypeMap, GenotypeMatrix> take_dominance_matrix() &&
    {
        return std::move(*dominance_matrix_);
    }
    bool has_dominance_matrix() const { return dominance_matrix_ != nullptr; }

    const std::vector<std::string>& fixed_effect_names() const;

   private:
    DataPipe() = default;

    template <typename Processor, typename TargetPtr>
    void load_genotype_impl(
        const std::string& suffix,
        TargetPtr& target,
        std::function<void(size_t, size_t)> progress_callback)
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
            auto result = pipe.template process<Processor>(
                config_.chunk_size, progress_callback);
            assign_to_target(std::move(result));
        }
        else
        {
            auto loader
                = gelex::GenotypeLoader(config_.bed_path, sample_manager_);
            auto result = loader.template process<Processor>(
                config_.chunk_size, progress_callback);
            assign_to_target(std::move(result));
        }
    }

    void convert_to_matrices();

    Config config_;

    std::unique_ptr<detail::PhenotypeLoader> phenotype_loader_;
    std::unique_ptr<detail::QcovarLoader> qcovar_loader_;
    std::unique_ptr<detail::DcovarLoader> dcovar_loader_;

    Eigen::VectorXd phenotype_;
    Eigen::MatrixXd fixed_effects_;

    std::shared_ptr<SampleManager> sample_manager_;

    std::unique_ptr<std::variant<GenotypeMap, GenotypeMatrix>> additive_matrix_;
    std::unique_ptr<std::variant<GenotypeMap, GenotypeMatrix>>
        dominance_matrix_;

    std::vector<std::string> fixed_effect_names_;
};

}  // namespace gelex
