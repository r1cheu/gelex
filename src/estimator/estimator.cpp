#include "gelex/estimator/estimator.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <armadillo>

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
        join_formula(model.random().random_indices(), model.random(), " + "),
        join_formula(model.random().genetic_indices(), model.random(), " + "),
        join_formula(model.random().gxe_indices(), model.random(), " + "));
    logger_->info(" \u25AA Samples:  {:d}", model.n_individuals());
    logger_->info("");

    logger_->info(fmt::format("{}", wine_red("[Term Summary]")));
    logger_->info(" \u25AA Fixed:  {}", fmt::join(model.fixed().names, ", "));
    if (model.random().has_random_effects())
    {
        logger_->info(
            " \u25AA Random:  {}",
            join_name(model.random().random_indices(), model.random(), ", "));
    }
    if (model.random().has_genetic_effects())
    {
        logger_->info(
            " \u25AA Genetic:  {}",
            join_name(model.random().genetic_indices(), model.random(), ", "));
    }
    if (model.random().has_gxe_effects())
    {
        logger_->info(
            " \u25AA GxE:  {}",
            join_name(model.random().gxe_indices(), model.random(), ", "));
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
        join_variance(model.random()),
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
            rebecca_purple(fmt::join(model.random().sigma(), ", ")),
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
    size_t iter{1};
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
            pink(fmt::join(model.random().sigma(), " ")),
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
    model.random().set_se(compute_se(model.random().hess_inv()));
    logger_->info(
        "───────────────────────────── Result "
        "───────────────────────────────");
    auto end = std::chrono::steady_clock::now();
    model.fixed().beta = compute_beta(model);
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
    for (size_t i = 0; i < model.n_fixed_effects(); ++i)
    {
        logger_->info(
            " \u25AA {}:  {:.6f} \u00B1 {:.4f}",
            model.fixed().levels[i],
            pink(model.fixed().beta.at(i)),
            pink(fixed_effse[i]));
    }
    logger_->info("");
}

void Estimator::report_variance_components(const GBLUP& model)
{
    logger_->info(fmt::format("{}", wine_red("[Variance Componests]")));
    report_variance("Random", model.random().random_indices(), model);
    report_variance("Genetic", model.random().genetic_indices(), model);
    report_variance("GxE", model.random().gxe_indices(), model);
    logger_->info(" \u25AA Residual:");
    logger_->info(
        "  - e:  {:6f} \u00B1 {:.4f}",
        pink(model.random().get("e")->sigma),
        pink(model.random().get("e")->se));
    logger_->info("");
}

void Estimator::report_heritability(const GBLUP& model)
{
    logger_->info(fmt::format("{}", wine_red("[Hertiability]")));
    auto [h2_se, sum_var] = compute_h2_se(model.random());
    size_t index{};
    for (auto genetic_index : model.random().genetic_indices())
    {
        logger_->info(
            " \u25AA {}:  {:.4f} \u00B1 {:.4f}",
            model.random()[genetic_index].name,
            pink(model.random()[genetic_index].sigma / sum_var),
            pink(h2_se[index]));
        ++index;
    }
    logger_->info("");
}

dvec Estimator::compute_beta(GBLUP& model)
{
    auto compute_beta_visitor = [this, &model](const auto& mat) -> dvec
    {
        return arma::inv_sympd(optimizer_->tx_vinv_x())
               * (mat.t() * optimizer_->v() * model.phenotype());
    };
    return std::visit(compute_beta_visitor, model.fixed().design_mat);
}

dmat Estimator::compute_u(GBLUP& model)
{
    dmat U{
        model.n_individuals(),
        model.random().sigma().n_elem,
        arma::fill::zeros};

    size_t idx{};
    for (const auto& effect : model.random())
    {
        std::visit(
            [&](const auto& cov)
            { U.unsafe_col(idx) = cov * optimizer_->proj_y() * effect.sigma; },
            effect.cov_mat);
    }
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
            model.random()[i].name,
            pink(model.random()[i].sigma),
            pink(model.random()[i].se));
    }
}

double Estimator::compute_aic(GBLUP& model)
{
    auto k
        = static_cast<double>(model.random().size() + model.n_fixed_effects());
    double aic = (-2 * optimizer_->loglike()) + (2 * k);
    return aic;
}

double Estimator::compute_bic(GBLUP& model)
{
    auto k
        = static_cast<double>(model.random().size() + model.n_fixed_effects());
    auto n = static_cast<double>(model.n_individuals());
    double bic = (-2 * optimizer_->loglike()) + (k * std::log(n));
    return bic;
}

dvec compute_se(const dmat& hess_inv)
{
    return arma::sqrt(arma::diagvec(-hess_inv));
}

std::pair<std::vector<double>, double> compute_h2_se(
    const RandomEffectManager& effects)
{
    auto n = effects.size();
    double sum_var = arma::sum(effects.sigma());
    double sum_sq = sum_var * sum_var;

    dvec grad = arma::zeros<dvec>(n);
    std::vector<double> h2_se;
    h2_se.reserve(effects.n_genetic_effects());
    for (auto i : effects.genetic_indices())
    {
        for (size_t j{}; j < n; ++j)
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
            std::sqrt(arma::as_scalar(grad.t() * -effects.hess_inv() * grad)));
    }
    return {h2_se, sum_var};
}

std::string join_formula(
    const std::vector<size_t>& indices,
    const RandomEffectManager& effects,
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
    const std::vector<size_t>& indices,
    const RandomEffectManager& effects,
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

std::string join_variance(const RandomEffectManager& effects)
{
    std::string result;
    for (const auto& effect : effects)
    {
        result += fmt::format(" V[{}] ", effect.name);
    }
    return result;
}

}  // namespace gelex
