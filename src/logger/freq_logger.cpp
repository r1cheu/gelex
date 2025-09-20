#include "../src/logger/freq_logger.h"

#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <armadillo>

#include "../src/logger/logger_utils.h"
#include "../src/model/freq/freq_effects.h"
#include "gelex/logger.h"
#include "gelex/model/freq/model.h"

namespace gelex
{
namespace detail
{

using arma::dvec;
EstimatorLogger::EstimatorLogger() : logger_{gelex::logging::get()} {}

void EstimatorLogger::set_verbose(bool verbose)
{
    if (!verbose)
    {
        logger_->set_level(spdlog::level::warn);
    }
}

void EstimatorLogger::log_model_info(
    const GBLUP& model,
    std::string_view optimizer_name,
    double tol,
    size_t max_iter)
{
    logger_->info(detail::title(" GBLUP MODEL ANALYSIS "));
    logger_->info(detail::subtitle("Model Specification"));
    logger_->info(detail::item("Model:  {}", model.formula()));
    logger_->info(detail::item("Samples:  {:d}", model.n_individuals()));
    logger_->info("");

    logger_->info(detail::subtitle("Optimizer Specification"));
    logger_->info(detail::item("Method:  {}", detail::cyan(optimizer_name)));
    logger_->info(detail::item("tolerance:  {:.2e}", tol));
    logger_->info(detail::item("Max Iterations:  {:d}", max_iter));
    logger_->info("");
}

void EstimatorLogger::log_em_initialization(
    double loglike,
    const TotalEffects& effects,
    double time_cost)
{
    logger_->info("Initializing with {} algorithm", detail::cyan("EM"));
    logger_->info(
        "Initial: logL={:.3f} | \u03C3\u00B2=[{:.3f}] ({:.3f}s)",
        loglike,
        detail::rebecca_purple(fmt::join(effects.values(), ", ")),
        time_cost);
}

void EstimatorLogger::log_iter_header(const GBLUP& model)
{
    logger_->info(detail::title(" REML ESTIMATION "));
    logger_->info(
        "{:>9} {:>9} {} {:>9}",
        "Iter.",
        "logL",
        join_variance(model.effects()),
        "duration");
}

void EstimatorLogger::log_iteration(
    size_t iter,
    double loglike,
    const TotalEffects& effects,
    double time_cost)
{
    logger_->info(
        "{:>9} {:>9.3f} "
        "{:>9.3f} {:>9.3f}s",
        iter,
        loglike,
        fmt::join(effects.values(), " "),
        time_cost);
}

void EstimatorLogger::log_results_header()
{
    logger_->info(detail::title(" RESULT "));
}

void EstimatorLogger::log_convergence_status(
    bool converged,
    size_t iter_count,
    size_t max_iter,
    double elapsed_time,
    double aic,
    double bic)
{
    logger_->info(detail::subtitle("Convergence"));

    if (converged)
    {
        logger_->info(
            " \u25AA Status:  {} ({} iterations in {:.3f}s)",
            detail::green("Success"),
            iter_count,
            elapsed_time);
    }
    else
    {
        logger_->warn(
            " \u25AA Status:  {} ({} iterations in {:.3f}s)",
            detail::red("Failed"),
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
    logger_->info(detail::subtitle("Fixed Effects"));
    for (size_t i = 0; i < model.n_fixed_effects(); ++i)
    {
        logger_->info(
            detail::item(
                "{}:  {}",
                model.fixed().levels[i],
                detail::with_std(model.fixed().coeff.at(i), fixed_se[i])));
    }
    logger_->info("");
}

void EstimatorLogger::log_variance_components(const GBLUP& model)
{
    logger_->info(detail::subtitle("Variance Components"));

    for (const auto& effect : model.effects())
    {
        logger_->info(
            detail::item(
                "{}:  {}",
                effect.name,
                detail::with_std(effect.sigma, effect.se)));
    }
    logger_->info("");
}

void EstimatorLogger::log_heritability(
    const GBLUP& model,
    const std::vector<double>& h2_se,
    double sum_var)
{
    logger_->info(detail::subtitle("Heritability"));
    size_t index{};
    for (const auto& effect : model.effects())
    {
        if (effect.type != effect_type::genetic)
        {
            continue;
        }

        logger_->info(
            detail::item(
                "{}:  {}",
                effect.name,
                detail::with_std(effect.sigma / sum_var, h2_se[index])));
        ++index;
    }
}

void EstimatorLogger::log_results_footer()
{
    logger_->info(detail::title(""));
}

std::string join_formula(
    const std::vector<size_t>& indices,
    const freq::RandomEffects& effects,
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
    const freq::RandomEffects& effects,
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

std::string join_variance(const TotalEffects& effects)
{
    std::string result;

    for (const auto& effect : effects)
    {
        result += fmt::format("{:>9}", fmt::format("V[{}]", effect.name));
    }
    return result;
}

}  // namespace detail
}  // namespace gelex
