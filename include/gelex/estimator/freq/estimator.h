#ifndef GELEX_ESTIMATOR_FREQ_ESTIMATOR_H_
#define GELEX_ESTIMATOR_FREQ_ESTIMATOR_H_

#include <cstddef>
#include <memory>

#include <spdlog/logger.h>

#include "gelex/optim/optimizer.h"
#include "gelex/optim/optimizer_state.h"

namespace gelex
{

class FreqModel;
class FreqState;

class Estimator
{
   public:
    explicit Estimator(size_t max_iter = 100, double tol = 1e-8);

    auto fit(
        const FreqModel& model,
        FreqState& state,
        bool em_init = true,
        bool verbose = true) -> Eigen::MatrixXd;

    auto is_converged() const -> bool { return converged_; }
    auto iter_count() const -> size_t { return iter_count_; }
    auto loglike() const -> double { return loglike_; }

   private:
    auto em_step(
        const FreqModel& model,
        FreqState& state,
        OptimizerState& opt_state) -> void;

    auto report_results(
        const FreqModel& model,
        const FreqState& state,
        const OptimizerState& opt_state,
        double elapsed) -> void;

    Optimizer optimizer_;
    size_t max_iter_{100};
    size_t iter_count_{};
    double tol_{1e-8};
    double loglike_{};
    bool converged_{false};

    std::shared_ptr<spdlog::logger> logger_;
};

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_FREQ_ESTIMATOR_H_
