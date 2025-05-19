#pragma once
#include <cstddef>
#include <memory>
#include <string_view>

#include <armadillo>

#include <fmt/ranges.h>
#include <spdlog/logger.h>
#include "gelex/model/effects.h"
#include "gelex/model/gblup.h"
#include "gelex/optim/base_optimizer.h"

namespace gelex
{

class Estimator
{
   public:
    Estimator(std::string_view optimizer, size_t max_iter, double tol);
    void set_optimizer(std::string_view optimizer, double tol);
    void fit(GBLUP& model, bool em_init = true, bool verbose = true);

   private:
    dvec compute_beta(GBLUP& model);
    dmat compute_u(GBLUP& model);
    std::unique_ptr<OptimizerBase> optimizer_;
    std::shared_ptr<spdlog::logger> logger_;

    void log_model_information(const GBLUP& model);
    void initialize_optimizer(GBLUP& model, bool em_init);
    void run_optimization_loop(GBLUP& model);

    void report_results(
        GBLUP& model,
        const std::chrono::steady_clock::time_point& start_time);
    void report_convergence_status(
        GBLUP& model,
        const std::chrono::steady_clock::time_point& start_time,
        const std::chrono::steady_clock::time_point& end_time);
    void report_fixed_effects(const GBLUP& model);
    void report_variance_components(const GBLUP& model);
    void report_heritability(const GBLUP& model);
    void report_variance(
        const std::string& category,
        const std::vector<size_t>& indices,
        const GBLUP& model);
    double compute_aic(GBLUP& model);
    double compute_bic(GBLUP& model);

    size_t iter_count_{};
    size_t max_iter_{};
    double tol_{};
    std::string optimizer_name_;
    bool converged_{};
};

template <typename eT>
auto blue_vec(const arma::Col<eT>& vec)
{
    return fmt::styled(fmt::join(vec, ", "), fmt::fg(fmt::color::blue_violet));
}

std::string join_formula(
    const std::vector<size_t>& indices,
    const RandomEffectManager& effects,
    std::string_view sep);

std::string join_name(
    const std::vector<size_t>& indices,
    const RandomEffectManager& effects,
    std::string_view sep);

std::string join_variance(const RandomEffectManager& effects);

dvec compute_se(const dmat& hess_inv);
std::pair<std::vector<double>, double> compute_h2_se(
    const RandomEffectManager& effects);

};  // namespace gelex
