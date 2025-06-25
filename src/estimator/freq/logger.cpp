#include "gelex/estimator/freq/logger.h"

#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <armadillo>

#include "gelex/utils/formatter.h"
#include "gelex/utils/utils.h"

namespace gelex
{

using arma::dvec;
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
    logger_->info(title(" GBLUP MODEL ANALYSIS "));
    logger_->info(subtitle("Model Specification"));
    logger_->info(item("Model:  {}", model.formula()));
    logger_->info(item("Samples:  {:d}", model.n_individuals()));
    logger_->info("");

    logger_->info(subtitle("Optimizer Specification"));
    logger_->info(item("Method:  {}", cyan(optimizer_name)));
    logger_->info(item("tolerance:  {:.2e}", tol));
    logger_->info(item("Max Iterations:  {:d}", max_iter));
    logger_->info(title(" REML ESTIMATION "));
    logger_->info(
        "{:>9} {:>9} {} {:>9}",
        "Iter.",
        "logL",
        join_variance(model.random()),
        "duration");
}

void EstimatorLogger::log_em_initialization(
    double loglike,
    const RandomEffectManager& effects,
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
    const RandomEffectManager& effects,
    double time_cost)
{
    logger_->info(
        "{:>9} {:>9.3f} "
        "{:>9.3f} {:>9.3f}s",
        iter,
        loglike,
        fmt::join(effects.sigma(), " "),
        time_cost);
}

void EstimatorLogger::log_results_header()
{
    logger_->info(title(" RESULT "));
}

void EstimatorLogger::log_convergence_status(
    bool converged,
    size_t iter_count,
    size_t max_iter,
    double elapsed_time,
    double aic,
    double bic)
{
    logger_->info(subtitle("Convergence"));

    if (converged)
    {
        logger_->info(
            " \u25AA Status:  {} ({} iterations in {:.3f}s)",
            green("Success"),
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
    logger_->info(" \u25AA AIC:  {:.3f}", aic);
    logger_->info(" \u25AA BIC:  {:.3f}", bic);
    logger_->info("");
}

void EstimatorLogger::log_fixed_effects(
    const GBLUP& model,
    const dvec& fixed_se)
{
    logger_->info(subtitle("Fixed Effects"));
    for (size_t i = 0; i < model.n_fixed_effects(); ++i)
    {
        logger_->info(item(
            "{}:  {}",
            model.fixed().levels[i],
            with_std(model.fixed().beta.at(i), fixed_se[i])));
    }
    logger_->info("");
}

void EstimatorLogger::log_variance_components(const GBLUP& model)
{
    logger_->info(subtitle("Variance Components"));
    log_variance_category("Random", model.random().random_indices(), model);
    log_variance_category("Genetic", model.random().genetic_indices(), model);
    log_variance_category("GxE", model.random().gxe_indices(), model);
    logger_->info(item("Residual:"));

    logger_->info(subitem(
        "e:  {}",
        with_std(model.random().get("e")->sigma, model.random().get("e")->se)));
    logger_->info("");
}

void EstimatorLogger::log_heritability(
    const GBLUP& model,
    const std::vector<double>& h2_se,
    double sum_var)
{
    logger_->info(subtitle("Heritability"));
    size_t index{};
    for (auto genetic_index : model.random().genetic_indices())
    {
        logger_->info(item(
            "{}:  {}",
            model.random()[genetic_index].name,
            with_std(
                model.random()[genetic_index].sigma / sum_var, h2_se[index])));
        ++index;
    }
}

void EstimatorLogger::log_results_footer()
{
    logger_->info(title(""));
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
        logger_->info(subitem(
            "{}:  {}",
            model.random()[i].name,
            with_std(model.random()[i].sigma, model.random()[i].se)));
    }
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
        result += fmt::format("{:>9}", fmt::format("V[{}]", effect.name));
    }
    return result;
}

}  // namespace gelex
