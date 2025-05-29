#pragma once
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <spdlog/logger.h>
#include "gelex/model/freq_effects.h"
#include "gelex/model/gblup.h"

namespace gelex
{

class EstimatorLogger
{
   public:
    EstimatorLogger();

    void set_verbose(bool verbose);
    void log_model_information(
        const GBLUP& model,
        std::string_view optimizer_name,
        double tol,
        size_t max_iter);
    void log_em_initialization(
        double loglike,
        const RandomEffectManager& effects,
        double time_cost);
    void log_iteration(
        size_t iter,
        double loglike,
        const RandomEffectManager& effects,
        double time_cost);
    void log_results_header();
    void log_convergence_status(
        bool converged,
        size_t iter_count,
        size_t max_iter,
        double elapsed_time,
        double aic,
        double bic);
    void log_fixed_effects(const GBLUP& model, const dvec& fixed_se);
    void log_variance_components(const GBLUP& model);
    void log_heritability(
        const GBLUP& model,
        const std::vector<double>& h2_se,
        double sum_var);
    void log_results_footer();

   private:
    std::shared_ptr<spdlog::logger> logger_;

    void log_variance_category(
        const std::string& category,
        const std::vector<size_t>& indices,
        const GBLUP& model);
};

// Helper functions for formatting
std::string join_formula(
    const std::vector<size_t>& indices,
    const RandomEffectManager& effects,
    std::string_view sep);

std::string join_name(
    const std::vector<size_t>& indices,
    const RandomEffectManager& effects,
    std::string_view sep);

std::string join_variance(const RandomEffectManager& effects);

}  // namespace gelex
