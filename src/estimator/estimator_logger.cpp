#include "gelex/estimator/estimator_logger.h"

#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <armadillo>

#include "gelex/utils.h"

namespace gelex
{

EstimatorLogger::EstimatorLogger() : logger_{Logger::logger()} {}

void EstimatorLogger::set_verbose(bool verbose)
{
    if (!verbose)
    {
        logger_->set_level(spdlog::level::warn);
    }
}

void EstimatorLogger::log_model_information(
    const GBLUP& model,
    std::string_view optimizer_name,
    double tol,
    size_t max_iter)
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
    logger_->info(" \u25AA Method:  {}", cyan(optimizer_name));
    logger_->info(" \u25AA tolerance:  {:.2e}", tol);
    logger_->info(" \u25AA Max Iterations:  {:d}", max_iter);
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

void EstimatorLogger::log_em_initialization(
    double loglike,
    const freq::RandomEffectManager& effects,
    double time_cost)
{
    logger_->info("Initializing with {} algorithm", cyan("EM"));
    logger_->info(
        "Initial: logL={:.3f} | \u03C3\u00B2=[{:.3f}] ({:.3f}s)",
        loglike,
        rebecca_purple(fmt::join(effects.sigma(), ", ")),
        time_cost);
}

void EstimatorLogger::log_iteration(
    size_t iter,
    double loglike,
    const freq::RandomEffectManager& effects,
    double time_cost)
{
    logger_->info(
        " {:^2d}{:>12.3f}  "
        "{:>7.3f} ({:.3f}s)",
        iter,
        pink(loglike),
        pink(fmt::join(effects.sigma(), " ")),
        time_cost);
}

void EstimatorLogger::log_results_header()
{
    logger_->info(
        "───────────────────────────── Result "
        "───────────────────────────────");
}

void EstimatorLogger::log_convergence_status(
    bool converged,
    size_t iter_count,
    size_t max_iter,
    double elapsed_time,
    double aic,
    double bic)
{
    logger_->info(fmt::format("{}", wine_red("[Convergence]")));

    if (converged)
    {
        logger_->info(
            " \u25AA Status:  {} ({} iterations in {:.3f}s)",
            fmt::format("{}", green("Success")),
            iter_count,
            elapsed_time);
    }
    else
    {
        logger_->warn(
            " \u25AA Status:  {} ({} iterations in {:.3f}s)",
            red("Failed"),
            max_iter,
            elapsed_time);
        logger_->warn(
            "Try to increase the max_iter or check the model specification.");
    }
    logger_->info(" \u25AA AIC:  {:.3f}", pink(aic));
    logger_->info(" \u25AA BIC:  {:.3f}", pink(bic));
    logger_->info("");
}

void EstimatorLogger::log_fixed_effects(
    const GBLUP& model,
    const dvec& fixed_se)
{
    logger_->info(fmt::format("{}", wine_red("[Fixed Effects]")));
    for (size_t i = 0; i < model.n_fixed_effects(); ++i)
    {
        logger_->info(
            " \u25AA {}:  {:.6f} \u00B1 {:.4f}",
            model.fixed().levels[i],
            pink(model.fixed().beta.at(i)),
            pink(fixed_se[i]));
    }
    logger_->info("");
}

void EstimatorLogger::log_variance_components(const GBLUP& model)
{
    logger_->info(fmt::format("{}", wine_red("[Variance Componests]")));
    log_variance_category("Random", model.random().random_indices(), model);
    log_variance_category("Genetic", model.random().genetic_indices(), model);
    log_variance_category("GxE", model.random().gxe_indices(), model);
    logger_->info(" \u25AA Residual:");
    logger_->info(
        "  - e:  {:6f} \u00B1 {:.4f}",
        pink(model.random().get("e")->sigma),
        pink(model.random().get("e")->se));
    logger_->info("");
}

void EstimatorLogger::log_heritability(
    const GBLUP& model,
    const std::vector<double>& h2_se,
    double sum_var)
{
    logger_->info(fmt::format("{}", wine_red("[Hertiability]")));
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

void EstimatorLogger::log_results_footer()
{
    logger_->info(
        "──────────────────────────────────────────"
        "──────────────────────────");
}

void EstimatorLogger::log_variance_category(
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

std::string join_formula(
    const std::vector<size_t>& indices,
    const freq::RandomEffectManager& effects,
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
    const freq::RandomEffectManager& effects,
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

std::string join_variance(const freq::RandomEffectManager& effects)
{
    std::string result;
    for (const auto& effect : effects)
    {
        result += fmt::format(" V[{}] ", effect.name);
    }
    return result;
}

}  // namespace gelex
