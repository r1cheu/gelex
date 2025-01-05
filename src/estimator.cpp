#include "chenx/estimator.h"

#include <cmath>

#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>

#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <armadillo>

#include "chenx/model/linear_mixed_model.h"
#include "chenx/optim/base_optimizer.h"
#include "chenx/optim/expectation_maximization.h"
#include "chenx/optim/second_order_optimizer.h"
#include "chenx/utils.h"
namespace chenx
{
Estimator::Estimator(std::string_view optimizer, size_t max_iter, double tol)
    : logger_{Logger::logger()}
{
    set_optimizer(optimizer, max_iter, tol);
}

void Estimator::set_optimizer(
    std::string_view optimizer,
    size_t max_iter,
    double tol)
{
    std::string opt_lower = ToLowercase(optimizer);
    // Map of optimizer names to factory functions
    static const std::unordered_map<
        std::string,
        std::function<std::unique_ptr<OptimizerBase>(size_t, double)>>
        optimizers
        = {{"nr", CreateOptimizer<NewtonRaphsonOptimizer>},
           {"newtonraphson", CreateOptimizer<NewtonRaphsonOptimizer>},
           {"fs", CreateOptimizer<FisherScoringOptimizer>},
           {"fisherscoring", CreateOptimizer<FisherScoringOptimizer>},
           {"ai", CreateOptimizer<AverageInformationOptimizer>},
           {"averageinformation", CreateOptimizer<AverageInformationOptimizer>},
           {"em", CreateOptimizer<ExpectationMaximizationOptimizer>},
           {"expectationmaximization",
            CreateOptimizer<ExpectationMaximizationOptimizer>}};

    auto it = optimizers.find(opt_lower);
    if (it != optimizers.end())
    {
        optimizer_ = it->second(max_iter, tol);
    }
    else
    {
        throw std::invalid_argument(
            "Unknown optimizer: " + std::string(optimizer) + ", NR(NewtonRaphson), FS(FisherScoring), AI(AverageInformation) are supported.");
    }
}

void Estimator::Fit(LinearMixedModel& model, bool em_init)
{
    logger_->info(
        "Starting fit with {:d} samples, variance componets: {} at {:%Y-%m-%d "
        "%H:%M:%S}",
        model.y().n_elem,
        model.rand_names(),
        fmt::localtime(std::time(nullptr)));
    logger_->info(
        "Optimizer: [{}](max_iter={:d}, tol={:.3e})",
        cyan(optimizer_->name()),
        optimizer_->max_iter(),
        optimizer_->tol());
    if (em_init)
    {
        double time_cost{};
        ExpectationMaximizationOptimizer em_optimizer{};
        logger_->info("Using [{}] to initialize.", cyan(em_optimizer.name()));

        // Start timing
        {
            Timer timer{time_cost};
            model.set_sigma(em_optimizer.Step(model));
        }

        logger_->info(
            "Initial loglikelihood={:.2f}, variance componets=[{:.4f}] {:.3f}s",
            model.ComputeLogLikelihood(),
            fmt::styled(
                fmt::join(model.sigma(), " "),
                fmt::fg(fmt::color::blue_violet)),
            time_cost);
    }
    bool converged{optimizer_->Optimize(model)};

    if (converged)
    {
        model.set_beta(ComputeBeta(model));
        logger_->info(
            "{} beta=[{:.6f}], variance componets=[{:.6f}]",
            gold("Converged!!!"),
            fmt::styled(
                fmt::join(model.beta(), " "), fmt::fg(fmt::color::blue_violet)),
            fmt::styled(
                fmt::join(model.sigma(), " "),
                fmt::fg(fmt::color::blue_violet)));
    }
    else
    {
        logger_->warn(
            "{}, try to increase max_iter({:d}), beta=[{:.6f}], "
            "variance componets=[{:.6f}]",
            red("Not converged!!!"),
            optimizer_->max_iter(),
            fmt::styled(
                fmt::join(model.beta(), " "), fmt::fg(fmt::color::blue_violet)),
            fmt::styled(
                fmt::join(model.sigma(), " "),
                fmt::fg(fmt::color::blue_violet)));
    }
}

dvec Estimator::ComputeBeta(LinearMixedModel& model)
{
    return arma::inv_sympd(model.tx_vinv_x())
           * (model.X().t() * model.v() * model.y());
}
}  // namespace chenx
