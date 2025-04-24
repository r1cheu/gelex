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
#include <variant>

#include "gelex/model/effects.h"
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
            + ", AI(Average Information) is supported.");
    }
}

void Estimator::fit(GBLUP& model, bool em_init, bool verbose)
{
    auto start = std::chrono::steady_clock::now();
    if (!verbose)
    {
        logger_->set_level(spdlog::level::warn);
    }

    logger_->info(
        "─────────────────────── GBLUP MODEL ANALYSIS "
        "───────────────────────");
    logger_->info("");
    logger_->info("MODEL SPECIFICATION:");
    logger_->info(
        " \u25AA Model:  {} + {}e",
        model.formula(),
        join_formula(model.effect().genetic_indices(), model.effect(), " + "));

    logger_->info(" \u25AA Samples:  {:d}", model.n_individuals());
    logger_->info(" \u25AA Common Effects:  {:d}", model.n_common_effects());
    logger_->info(" \u25AA Group Effects:  {:d}", model.n_group_effects());
    logger_->info(
        " \u25AA Genetic Effects:  {}",
        join_name(model.effect().genetic_indices(), model.effect(), ", "));
    logger_->info("");
    logger_->info("OPTIMIZATION SPECIFICATION:");
    logger_->info(" \u25AA Method:  {}", cyan(optimizer_name_));
    logger_->info(" \u25AA tolerance:  {:.2e}", tol_);
    logger_->info(" \u25AA Max Iterations:  {:d}", max_iter_);
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
            " {:<2d}  logL={:.3f} | "
            "\u03C3\u00B2=[{:.3f}] ({:.3f}s)",
            iter,
            pink(optimizer_->loglike()),
            pink(fmt::join(model.sigma(), ", ")),
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
            " \u25AA Status:  {} ({} iterations in {:.3f}s)",
            gold("Success"),
            iter,
            std::chrono::duration<double>(end - start).count());
        logger_->info("");
    }
    else
    {
        logger_->warn(
            " \u25AA Status:  {} ({} iterations in {:.3f}s)",
            red("Failed"),
            max_iter_,
            std::chrono::duration<double>(end - start).count());
        logger_->warn(
            "Try to increase the max_iter or check the model specification.");
        logger_->info("");
    }

    logger_->info("COMMON EFFECTS: ");
    logger_->info(
        " \u25AA \u03B2:  [{:.6f}]", pink(fmt::join(model.beta(), ", ")));
    logger_->info("");

    logger_->info("VARIANCE COMPONESTS: ");
    report_variance("Group", model.effect().group_indices(), model);
    report_variance("Genetic", model.effect().genetic_indices(), model);
    report_variance("GxE", model.effect().gxe_indices(), model);
    logger_->info(" \u25AA Residual Variance:");
    logger_->info("  \u2022 e:  {:6f}", pink(model.sigma().back()));
    logger_->info("");

    logger_->info("HERTIABILITY:");
    const double sum_var = arma::sum(model.sigma());
    for (auto i : model.effect().genetic_indices())
    {
        logger_->info(
            " \u25AA {}:  {:.3f}",
            model.effect()[i].name,
            pink(model.effect()[i].sigma / sum_var));
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
    for (const Effect& effect : model.effect())
    {
        std::visit(
            [&](const auto& cov)
            { U.unsafe_col(idx) = cov * optimizer_->proj_y() * effect.sigma; },
            effect.cov_mat);
    }
    U.unsafe_col(idx) = optimizer_->proj_y() * model.sigma().back();
    return U;
}

void Estimator::report_variance(
    const std::string& category,
    const std::vector<size_t>& indices,
    const GBLUP& model)
{
    if (indices.empty())
    {
        return;
    }

    logger_->info(" \u25AA {} Variance:", category);
    for (auto i : indices)
    {
        logger_->info(
            "  \u2022 {}:  {:.6f}",
            model.effect()[i].name,
            pink(model.effect()[i].sigma));
    }
}

std::string join_formula(
    const std::vector<uint64_t>& indices,
    const GroupEffectManager& effects,
    std::string_view sep)
{
    std::string result;
    switch (effects[indices[0]].type)
    {
        case gelex::effect_type::group:
        {
            for (auto i : indices)
            {
                result += fmt::format(
                    "{}({}){}", effects[i].name, green("R[E]"), sep);
            }
            break;
        }
        case gelex::effect_type::genetic:
        {
            for (auto i : indices)
            {
                result += fmt::format(
                    "{}({}){}", effects[i].name, cyan("R[G]"), sep);
            }
            break;
        }
        case gelex::effect_type::gxe:
        {
            for (auto i : indices)
            {
                result += fmt::format(
                    "{}({}){}", effects[i].name, rebecca_purple("R[GxE]"), sep);
            }
            break;
        }
    };
    return result;
}

std::string join_name(
    const std::vector<uint64_t>& indices,
    const GroupEffectManager& effects,
    std::string_view sep)
{
    std::string result;
    for (size_t i = 0; i < indices.size(); ++i)
    {
        result += effects[indices[i]].name;
        if (i != indices.size() - 1)
        {
            result += sep;
        }
    }
    return result;
}
}  // namespace gelex
