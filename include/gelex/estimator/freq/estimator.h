#pragma once
#include <cstddef>
#include <string_view>

#include <fmt/ranges.h>
#include <armadillo>

#include "../src/logger/freq_logger.h"
#include "gelex/model/freq/model.h"
#include "gelex/optim/optimizer.h"

namespace gelex
{

class Estimator
{
   public:
    Estimator(std::string_view optimizer, size_t max_iter, double tol);
    void set_optimizer(std::string_view optimizer, double tol);
    void fit(GBLUP& model, bool em_init = true, bool verbose = true);

   private:
    void em_step(GBLUP& model, bool em_init);

    void compute_beta(GBLUP& model);
    void compute_u(GBLUP& model);
    void report_results(
        GBLUP& model,
        const std::chrono::steady_clock::time_point& start_time);

    double compute_aic(const GBLUP& model) const;
    double compute_bic(const GBLUP& model) const;
    void compute_var_se(GBLUP& model, const Optimizer& optim) const;
    std::pair<std::vector<double>, double> compute_h2_se(
        const GBLUP& model) const;

    Optimizer optimizer_;
    detail::EstimatorLogger logger_;

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

};  // namespace gelex
