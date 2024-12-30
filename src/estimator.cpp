#include "chenx/estimator.h"
#include <spdlog/spdlog.h>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include "armadillo"

#include "chenx/logger.h"
#include "chenx/optim/expectation_maximization.h"
#include "chenx/optim/second_order_optimizer.h"
#include "chenx/utils.h"
namespace chenx
{
Estimator::Estimator(std::string_view optimizer, size_t max_iter, double tol)
    : logger_{spdlog::get(logger_name)}
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
    logger_->info("Optimize using [{}]", optimizer_->name());
}

void Estimator::Fit(LinearMixedModel& model, bool em_init)
{
    if (em_init)
    {
        logger_->info("Using [ExpectationMaximization] to initialize");
        ExpectationMaximizationOptimizer em_optimizer{};
        model.set_sigma(em_optimizer.Step(model));
    }
    bool converged{optimizer_->Optimize(model)};

    if (converged)
    {
        model.set_beta(ComputeBeta(model));
    }
}

dvec Estimator::ComputeBeta(LinearMixedModel& model)
{
    return arma::inv_sympd(model.tx_vinv_x())
           * (model.X().t() * model.v() * model.y());
}
}  // namespace chenx
