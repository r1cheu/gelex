#include "gelex/estimator/estimator.h"

#include <cmath>

#include <cstdint>
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
    auto start = std::chrono::steady_clock::now();
    if (!verbose)
    {
        logger_->set_level(spdlog::level::warn);
    }

    model.set_model();

    logger_->info(
        "─────────────────────── GBLUP MODEL ANALYSIS "
        "───────────────────────");
    logger_->info("");
    logger_->info("MODEL SPECIFICATION:");
    logger_->info(
        "\u2022 Model:     {} + {} + e",
        model.formula(),
        fmt::join(model.genetic_names(), " + "));

    logger_->info("\u2022 Samples:     {:d}", model.n_individuals());
    logger_->info("\u2022 Common Effects: {:d}", model.n_common_effects());
    logger_->info("\u2022 Group Effects: {:d}", model.n_group_effects());
    logger_->info(
        "\u2022 Genetic Effects: {}", fmt::join(model.genetic_names(), ", "));
    logger_->info("");
    logger_->info("OPTIMIZATION SPECIFICATION:");
    logger_->info("\u2022 Method:    {}", cyan(optimizer_name_));
    logger_->info("\u2022 tolerance:    {:.2e}", tol_);
    logger_->info("\u2022 Max Iterations:   {:d}", max_iter_);
    logger_->info("");
    logger_->info(
        "────────────────────────── REML ITERATIONS "
        "─────────────────────────");
    double time_cost{};
    if (em_init)
    {
        ExpectationMaximizationOptimizer em_optimizer{tol_};
        logger_->info("Initializing with {} algorithm", cyan("EM"));
        em_optimizer.init(model);

        {
            Timer timer{time_cost};
            em_optimizer.step(model);
        }

        logger_->info(
            "Initial: logL={:.3f} | \u03C3\u00B2=[{:.3f}] ({:.3f}s)",
            em_optimizer.loglike(),
            rebecca_purple(fmt::join(model.sigma(), ", ")),
            time_cost);
        optimizer_ = std::make_unique<AverageInformationOptimizer>(
            std::move(static_cast<OptimizerBase&>(em_optimizer)));
    }
    else
    {
        optimizer_ = std::make_unique<AverageInformationOptimizer>(tol_);
        optimizer_->init(model);
    }
    uint64_t iter{1};
    for (; iter <= max_iter_; ++iter)
    {
        {
            Timer timer{time_cost};
            optimizer_->step(model);
        }
        logger_->info(
            "{:<2d}  logL={:.3f} | "
            "\u03C3\u00B2=[{:.3f}] ({:.3f}s)",
            iter,
            optimizer_->loglike(),
            rebecca_purple(fmt::join(model.sigma(), ", ")),
            time_cost);

        if (optimizer_->converged())
        {
            converged_ = true;
            break;
        }
    }
    logger_->info(
        "───────────────────────────── Result "
        "───────────────────────────────");
    auto end = std::chrono::steady_clock::now();
    model.set_beta(compute_beta(model));
    logger_->info("CONVERGENCE:");

    if (converged_)
    {
        logger_->info(
            "\u2022 Status:      {} ({} iterations in {:.3f}s)",
            gold("Success"),
            iter,
            std::chrono::duration<double>(end - start).count());
        logger_->info("");
    }
    else
    {
        logger_->warn(
            "\u2022 Status:      {} ({} iterations in {:.3f}s)",
            red("Failed"),
            max_iter_,
            std::chrono::duration<double>(end - start).count());
        logger_->warn(
            "Try to increase the max_iter or check the model specification.");
        logger_->info("");
    }

    logger_->info("COMMON EFFECTS: ");
    logger_->info(
        "\u2022 \u03B2:    [{:.6f}]",
        rebecca_purple(fmt::join(model.beta(), ", ")));
    logger_->info("");
    logger_->info("VARIANCE COMPONESTS: ");
    for (size_t i = 0; i < model.n_group_effects(); ++i)
    {
        logger_->info(
            "\u2022 {}:     {:.6f}",
            model.group_names()[i],
            rebecca_purple(model.sigma().at(i)));
    }
    for (size_t i = 0; i < model.n_genetic_effects(); ++i)
    {
        logger_->info(
            "\u2022 {}:     {:.6f}",
            model.genetic_names()[i],
            rebecca_purple(model.sigma().at(i + model.n_group_effects())));
    }
    logger_->info("\u2022 e:     {:.6f}", rebecca_purple(model.sigma().back()));

    logger_->info("");
    logger_->info("HERTIABILITY:");

    dvec h2 = compute_h2(model);
    for (size_t i{0}; i < model.n_genetic_effects(); ++i)
    {
        logger_->info(
            "\u2022 h\u00B2_{}:     {:.3f}",
            model.genetic_names()[i],
            rebecca_purple(h2.at(i)));
    }
    logger_->info(
        "──────────────────────────────────────────"
        "──────────────────────────");
    compute_u(model);
}

dvec Estimator::compute_beta(GBLUP& model)
{
    return arma::inv_sympd(optimizer_->tx_vinv_x())
           * (model.design_mat_beta().t() * optimizer_->v()
              * model.phenotype());
}

dmat Estimator::compute_u(GBLUP& model)
{
    dmat U{model.n_individuals(), model.sigma().n_elem, arma::fill::zeros};

    uint64_t idx{};
    for (const sp_dmat& mat : model.group_cov_mats())
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
