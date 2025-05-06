#include "gelex/estimator/estimator.h"

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

    log_model_information(model);
    initialize_optimizer(model, em_init);
    run_optimization_loop(model);
    report_results(model, start);
}

void Estimator::log_model_information(const GBLUP& model)
{
    logger_->info(
        "─────────────────────── GBLUP MODEL ANALYSIS "
        "───────────────────────");
    logger_->info("");
    logger_->info(fmt::format("{}", wine_red("[Model Specification]")));
    logger_->info(
        " \u25AA Model:  {}{}{}{}e",
        model.formula(),
        join_formula(model.effect().group_indices(), model.effect(), " + "),
        join_formula(model.effect().genetic_indices(), model.effect(), " + "),
        join_formula(model.effect().gxe_indices(), model.effect(), " + "));
    logger_->info(" \u25AA Samples:  {:d}", model.n_individuals());
    logger_->info("");

    logger_->info(fmt::format("{}", wine_red("[Term Summary]")));
    logger_->info(
        " \u25AA Fixed:  {}", fmt::join(model.fixed_effect_names(), ", "));
    if (model.effect().has_group_effects())
    {
        logger_->info(
            " \u25AA Random:  {}",
            join_name(model.effect().group_indices(), model.effect(), ", "));
    }
    if (model.effect().has_genetic_effects())
    {
        logger_->info(
            " \u25AA Genetic:  {}",
            join_name(model.effect().genetic_indices(), model.effect(), ", "));
    }
    if (model.effect().has_gxe_effects())
    {
        logger_->info(
            " \u25AA GxE:  {}",
            join_name(model.effect().gxe_indices(), model.effect(), ", "));
    }
    logger_->info("");

    logger_->info(fmt::format("{}", wine_red("[Optimizer Specification]")));
    logger_->info(" \u25AA Method:  {}", cyan(optimizer_name_));
    logger_->info(" \u25AA tolerance:  {:.2e}", tol_);
    logger_->info(" \u25AA Max Iterations:  {:d}", max_iter_);
    logger_->info("");
    logger_->info(
        "────────────────────────── REML ITERATIONS "
        "─────────────────────────");
    logger_->info(
        "{:>5}{:>8}  {}  {:<10}",
        "Iter.",
        "logL",
        join_variance(model.effect()),
        "duration");
}

void Estimator::initialize_optimizer(GBLUP& model, bool em_init)
{
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
}

void Estimator::run_optimization_loop(GBLUP& model)
{
    uint64_t iter{1};
    double time_cost{};

    for (; iter <= max_iter_; ++iter)
    {
        {
            Timer timer{time_cost};
            optimizer_->step(model);
        }
        logger_->info(
            " {:^2d}{:>12.3f}  "
            "{:>7.3f} ({:.3f}s)",
            iter,
            pink(optimizer_->loglike()),
            pink(fmt::join(model.sigma(), " ")),
            time_cost);

        if (optimizer_->converged())
        {
            converged_ = true;
            break;
        }
    }
    iter_count_ = iter;
}

void Estimator::report_results(
    GBLUP& model,
    const std::chrono::steady_clock::time_point& start_time)
{
    logger_->info(
        "───────────────────────────── Result "
        "───────────────────────────────");
    auto end = std::chrono::steady_clock::now();
    model.set_beta(compute_beta(model));
    logger_->info(fmt::format("{}", wine_red("[Convergence]")));

    report_convergence_status(model, start_time, end);
    report_fixed_effects(model);
    report_variance_components(model);
    report_heritability(model);

    logger_->info(
        "──────────────────────────────────────────"
        "──────────────────────────");
    compute_u(model);
}

void Estimator::report_convergence_status(
    GBLUP& model,
    const std::chrono::steady_clock::time_point& start_time,
    const std::chrono::steady_clock::time_point& end_time)
{
    double elapsed_time
        = std::chrono::duration<double>(end_time - start_time).count();

    if (converged_)
    {
        logger_->info(
            " \u25AA Status:  {} ({} iterations in {:.3f}s)",
            fmt::format("{}", green("Success")),
            iter_count_,
            elapsed_time);
    }
    else
    {
        logger_->warn(
            " \u25AA Status:  {} ({} iterations in {:.3f}s)",
            red("Failed"),
            max_iter_,
            elapsed_time);
        logger_->warn(
            "Try to increase the max_iter or check the model specification.");
    }
    logger_->info(" \u25AA AIC:  {:.3f}", pink(compute_aic(model)));
    logger_->info(" \u25AA BIC:  {:.3f}", pink(compute_bic(model)));
    logger_->info("");
}

void Estimator::report_fixed_effects(const GBLUP& model)
{
    logger_->info(fmt::format("{}", wine_red("[Fixed Effects]")));
    dvec fixed_effse
        = arma::diagvec(arma::sqrt(arma::inv_sympd(optimizer_->tx_vinv_x())));
    for (size_t i = 0; i < model.n_common_effects(); ++i)
    {
        logger_->info(
            " \u25AA {}:  {:.6f} \u00B1 {:.4f}",
            model.fixed_effect_levles()[i],
            pink(model.beta()[i]),
            pink(fixed_effse[i]));
    }
    logger_->info("");
}

void Estimator::report_variance_components(const GBLUP& model)
{
    logger_->info(fmt::format("{}", wine_red("[Variance Componests]")));
    report_variance("Random", model.effect().group_indices(), model);
    report_variance("Genetic", model.effect().genetic_indices(), model);
    report_variance("GxE", model.effect().gxe_indices(), model);
    logger_->info(" \u25AA Residual:");
    logger_->info(
        "  - e:  {:6f} \u00B1 {:.4f}",
        pink(model.sigma().back()),
        pink(model.effect().residual_se()));
    logger_->info("");
}

void Estimator::report_heritability(const GBLUP& model)
{
    logger_->info(fmt::format("{}", wine_red("[Hertiability]")));

    const double sum_var = arma::sum(model.sigma());
    std::vector<double> h2_se = compute_h2_se(sum_var, model.effect());
    uint64_t index{};
    for (auto genetic_index : model.effect().genetic_indices())
    {
        logger_->info(
            " \u25AA {}:  {:.4f} \u00B1 {:.4f}",
            model.effect()[genetic_index].name,
            pink(model.effect()[genetic_index].sigma / sum_var),
            pink(h2_se[index]));
        ++index;
    }
    logger_->info("");
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

    logger_->info(" \u25AA {}:", category);
    for (auto i : indices)
    {
        logger_->info(
            "  - {}:  {:.6f} \u00B1 {:.4f}",
            model.effect()[i].name,
            pink(model.effect()[i].sigma),
            pink(model.effect()[i].se));
    }
}

double Estimator::compute_aic(GBLUP& model)
{
    auto k
        = static_cast<double>(model.effect().size() + model.n_common_effects());
    double aic = (-2 * optimizer_->loglike()) + (2 * k);
    return aic;
}

double Estimator::compute_bic(GBLUP& model)
{
    auto k
        = static_cast<double>(model.effect().size() + model.n_common_effects());
    auto n = static_cast<double>(model.n_individuals());
    double bic = (-2 * optimizer_->loglike()) + (k * std::log(n));
    return bic;
}

std::vector<double> compute_h2_se(double sum_var, const EffectManager& effects)
{
    auto n = effects.size();
    double sum_sq = sum_var * sum_var;

    dvec grad = arma::zeros<dvec>(n);
    std::vector<double> h2_se;
    h2_se.reserve(effects.n_genetic_effects());
    for (auto i : effects.genetic_indices())
    {
        for (uint64_t j{}; j < n; ++j)
        {
            if (i == j)
            {
                grad[j] = (sum_var - effects[i].sigma) / sum_sq;
            }
            else
            {
                grad[j] = -effects[i].sigma / sum_sq;
            }
        }
        h2_se.emplace_back(
            std::sqrt(
                arma::as_scalar(grad.t() * -effects.effect_cov() * grad)));
    }
    return h2_se;
}

std::string join_formula(
    const std::vector<uint64_t>& indices,
    const EffectManager& effects,
    std::string_view sep)
{
    if (indices.empty())
    {
        return "";
    }

    std::string result;
    for (auto i : indices)
    {
        result += fmt::format("{}{}", effects[i].name, sep);
    }
    return result;
}

std::string join_name(
    const std::vector<uint64_t>& indices,
    const EffectManager& effects,
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

std::string join_variance(const EffectManager& effects)
{
    std::string result;
    for (const auto& effect : effects)
    {
        result += fmt::format(" V[{}] ", effect.name);
    }
    result += " V[e]";
    return result;
}

}  // namespace gelex
