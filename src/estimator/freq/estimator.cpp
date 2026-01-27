#include "gelex/estimator/freq/estimator.h"

#include <Eigen/Core>
#include <chrono>

#include "../src/logger/reml_logger.h"
#include "../src/logger/reml_logger_base.h"
#include "../src/utils/utils.h"
#include "gelex/estimator/freq/effect_solver.h"
#include "gelex/estimator/freq/statistics.h"
#include "gelex/model/freq/model.h"
#include "gelex/optim/policy.h"
#include "gelex/optim/variance_calculator.h"

namespace gelex
{

Estimator::Estimator(size_t max_iter, double tol)
    : optimizer_(tol),
      max_iter_(max_iter),
      logger_(std::make_unique<detail::RemlLogger>())
{
}

Estimator::Estimator(
    size_t max_iter,
    double tol,
    std::unique_ptr<detail::RemlLoggerBase> logger)
    : optimizer_(tol), max_iter_(max_iter), logger_(std::move(logger))
{
}

Estimator::~Estimator() = default;

auto Estimator::fit(
    const FreqModel& model,
    FreqState& state,
    bool em_init,
    bool verbose) -> Eigen::MatrixXd
{
    auto start = std::chrono::steady_clock::now();

    logger_->set_verbose(verbose);

    OptimizerState opt_state(model);

    // EM initialization
    if (em_init)
    {
        em_step(model, state, opt_state);
    }

    logger_->log_iter_header(state);

    // AI iterations
    for (size_t iter = 1; iter <= max_iter_; ++iter)
    {
        double time_cost{};
        {
            Timer timer{time_cost};
            optimizer_.step<AIPolicy>(model, state, opt_state);
        }

        loglike_ = variance_calculator::compute_loglike(model, opt_state);

        logger_->log_iteration(iter, loglike_, state, time_cost);

        if (optimizer_.is_converged())
        {
            converged_ = true;
            iter_count_ = iter;
            break;
        }
    }
    logger_->log_iter_footer();

    if (!converged_)
    {
        iter_count_ = max_iter_;
    }

    // compute final results
    effect_solver::compute_fixed_effects(model, state, opt_state);
    effect_solver::compute_random_effects(model, state, opt_state);
    statistics::compute_variance_se(state, opt_state);
    statistics::compute_heritability(state, opt_state);

    auto elapsed = std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - start)
                       .count();
    logger_->log_results(
        model, state, loglike_, converged_, iter_count_, max_iter_, elapsed);
    return std::move(opt_state.v);
}

auto Estimator::em_step(
    const FreqModel& model,
    FreqState& state,
    OptimizerState& opt_state) -> void
{
    double time_cost{};
    {
        Timer timer{time_cost};
        optimizer_.step<EMPolicy>(model, state, opt_state);
    }

    double loglike = variance_calculator::compute_loglike(model, opt_state);
    logger_->log_em_init(state, loglike);
}

}  // namespace gelex
