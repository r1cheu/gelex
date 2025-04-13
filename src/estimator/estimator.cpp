#include "gelex/estimator/estimator.h"

#include <cmath>

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <armadillo>

#include "gelex/model/gblup.h"
#include "gelex/optim/base_optimizer.h"
#include "gelex/optim/optimizers.h"
#include "gelex/utils.h"
namespace gelex
{
Estimator::Estimator(std::string_view optimizer, size_t max_iter, double tol)
    : logger_{Logger::logger()}, max_iter_{max_iter}
{
    set_optimizer(optimizer, tol);
}

void Estimator::set_optimizer(std::string_view optimizer, double tol)
{
    std::string opt_lower = ToLowercase(optimizer);
    // Map of optimizer names to factory functions
    if (opt_lower == "ai")
    {
        optimizer_name_ = "AI";
        tol_ = tol;
    }
    else
    {
        throw std::invalid_argument(
            "Unknown optimizer: " + std::string(optimizer)
            + "AI(Average Information) is supported.");
    }
}

void Estimator::fit(GBLUP& model, bool em_init, bool verbose)
{
    if (!verbose)
    {
        logger_->set_level(spdlog::level::warn);
    }

    model.set_model();

    logger_->info(
        "Start fitting with {:d} samples, {:d} common effects, variance "
        "componets: {} at {:%Y-%m-%d "
        "%H:%M:%S}",
        model.n_individuals(),
        model.n_common_effects(),
        fmt::join(model.sigma_names(), ", "),
        fmt::localtime(std::time(nullptr)));

    logger_->info(
        "Optimizer: [{}](tol={:.3e}), max_iter = {:d}",
        cyan(optimizer_name_),
        tol_,
        max_iter_);

    double time_cost{};
    if (em_init)
    {
        ExpectationMaximizationOptimizer em_optimizer{tol_};
        logger_->info("Using [{}] to initialize.", cyan("EM"));
        em_optimizer.init(model);

        // Start timing
        {
            Timer timer{time_cost};
            em_optimizer.step(model);
        }

        logger_->info(
            "Initial loglikelihood={:.2f}, variance componets=[{:.4f}] {:.3f}s",
            em_optimizer.loglike(),
            fmt::styled(
                fmt::join(model.sigma(), " "),
                fmt::fg(fmt::color::blue_violet)),
            time_cost);
        optimizer_ = std::make_unique<AverageInformationOptimizer>(
            std::move(static_cast<OptimizerBase&>(em_optimizer)));
    }
    else
    {
        optimizer_ = std::make_unique<AverageInformationOptimizer>(tol_);
        optimizer_->init(model);
    }

    for (uint64_t i{}; i < max_iter_; i++)
    {
        {
            Timer timer{time_cost};
            optimizer_->step(model);
            logger_->info(
                "Iter {:d}: logL={:.2f}, varcomp=[{:.4f}] {:.3f}s",
                i,
                optimizer_->loglike(),
                fmt::styled(
                    fmt::join(model.sigma(), " "),
                    fmt::fg(fmt::color::blue_violet)),
                time_cost);
        }

        if (optimizer_->converged())
        {
            converged_ = true;
            break;
        }
    }

    if (converged_)
    {
        model.set_beta(compute_beta(model));
        logger_->info(
            "{} beta=[{:.6f}], variance componets=[{:.6f}]",
            gold("Converged!!!"),
            fmt::styled(
                fmt::join(model.beta(), " "), fmt::fg(fmt::color::blue_violet)),
            fmt::styled(
                fmt::join(model.sigma(), " "),
                fmt::fg(fmt::color::blue_violet)));
        compute_u(model);
    }
    else
    {
        logger_->warn(
            "{}, try to increase max_iter({:d}), beta=[{:.6f}], "
            "variance componets=[{:.6f}]",
            red("Not converged!!!"),
            max_iter_,
            fmt::styled(
                fmt::join(model.beta(), " "), fmt::fg(fmt::color::blue_violet)),
            fmt::styled(
                fmt::join(model.sigma(), " "),
                fmt::fg(fmt::color::blue_violet)));
    }
}

dvec Estimator::compute_beta(GBLUP& model)
{
    return arma::inv_sympd(optimizer_->tx_vinv_x())
           * (model.design_mat_beta().t() * optimizer_->v()
              * model.phenotype());
}

dmat Estimator::compute_u(GBLUP& model)
{
    dmat U{model.n_individuals(), model.n_random_effects(), arma::fill::zeros};

    uint64_t idx;
    for (const sp_dmat& mat : model.env_cov_mats())
    {
        U.unsafe_col(idx) += mat * optimizer_->proj_y() * model.sigma().at(idx);
        ++idx;
    }
    for (const dmat& mat : model.genetic_cov_mats())
    {
        U.unsafe_col(idx) += mat * optimizer_->proj_y() * model.sigma().at(idx);
        ++idx;
    }
    U.unsafe_col(idx) = optimizer_->proj_y() * model.sigma().back();
    return U;
}
}  // namespace gelex
