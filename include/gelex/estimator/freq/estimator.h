#pragma once
#include <cstddef>
#include <memory>
#include <string_view>

#include <armadillo>

#include <fmt/ranges.h>
#include "gelex/estimator/freq/logger.h"
#include "gelex/model/freq/effects.h"
#include "gelex/model/freq/model.h"
#include "gelex/optim/base.h"

namespace gelex
{

class Estimator
{
   public:
    Estimator(std::string_view optimizer, size_t max_iter, double tol);
    void set_optimizer(std::string_view optimizer, double tol);
    void fit(GBLUP& model, bool em_init = true, bool verbose = true);

   private:
    arma::dvec compute_beta(GBLUP& model);
    arma::dmat compute_u(GBLUP& model);
    std::unique_ptr<OptimizerBase> optimizer_;
    EstimatorLogger logger_;

    void initialize_optimizer(GBLUP& model, bool em_init);
    void run_optimization_loop(GBLUP& model);

    void report_results(
        GBLUP& model,
        const std::chrono::steady_clock::time_point& start_time);

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

arma::dvec compute_se(const arma::dmat& hess_inv);
std::pair<std::vector<double>, double> compute_h2_se(
    const RandomEffectManager& effects);

};  // namespace gelex
