#include "../src/logger/reml_logger.h"

#include <fmt/format.h>
#include <fmt/ranges.h>

#include "../utils/formatter.h"
#include "gelex/estimator/freq/statistics.h"
#include "gelex/logger.h"
#include "gelex/model/freq/model.h"

namespace gelex
{
namespace detail
{

void RemlLogger::set_verbose(bool verbose)
{
    RemlLoggerBase::set_verbose(verbose);
}

void RemlLogger::log_em_init(const FreqState& state, double loglike)
{
    logger_->info(progress_mark("Initializing (EM)..."));

    std::string var_values;
    for (const auto& g : state.genetic())
    {
        var_values += fmt::format("{:.2f}, ", g.variance);
    }
    for (const auto& r : state.random())
    {
        var_values += fmt::format("{:.2f}, ", r.variance);
    }
    var_values += fmt::format("{:.2f}", state.residual().variance);

    logger_->info(
        "    LogL: {:.2f} | Init Vg: [{}]",
        loglike,
        rebecca_purple(var_values));
}

void RemlLogger::log_iter_header(const FreqState& state)
{
    logger_->info("");
    std::string var_header;
    for (const auto& g : state.genetic())
    {
        var_header += fmt::format("{:>12}", fmt::format("V({})", g.type));
    }
    for (const auto& r : state.random())
    {
        var_header += fmt::format("{:>12}", fmt::format("V({})", r.name));
    }
    var_header += fmt::format("{:>12}", "V(e)");

    logger_->info(
        "  {:<4} {:>12} {} {:>10}", "Iter", "LogL", var_header, "Time");
    logger_->info(gelex::table_separator(67));
}

void RemlLogger::log_iteration(
    size_t iter,
    double loglike,
    const FreqState& state,
    double time_cost)
{
    std::string var_values;
    for (const auto& g : state.genetic())
    {
        var_values += fmt::format("{:>12.2f}", g.variance);
    }
    for (const auto& r : state.random())
    {
        var_values += fmt::format("{:>12.2f}", r.variance);
    }
    var_values += fmt::format("{:>12.2f}", state.residual().variance);
    logger_->info(
        "  {:<4} {:>12.2f} {} {:>9.2f}s", iter, loglike, var_values, time_cost);
}

void RemlLogger::log_iter_footer()
{
    logger_->info(gelex::table_separator(67));
}

void RemlLogger::log_results(
    const FreqModel& model,
    const FreqState& state,
    double loglike,
    bool converged,
    size_t iter_count,
    size_t max_iter,
    double elapsed)
{
    logger_->info("");
    logger_->info(
        fmt::format(
            fmt::fg(fmt::color::light_cyan),
            "── REML Results {}",
            separator(70 - 16)));

    log_convergence(converged, iter_count, max_iter, elapsed);
    log_model_fit(model, loglike);
    log_fixed_effects(model, state);
    log_variance_components(state);

    logger_->info(
        fmt::format(fmt::fg(fmt::color::light_cyan), "{}", separator(70)));
}

void RemlLogger::log_convergence(
    bool converged,
    size_t iter_count,
    size_t max_iter,
    double elapsed)
{
    if (converged)
    {
        logger_->info(success(
            "Converged successfully in {} iterations ({:.2f}s)",
            iter_count,
            elapsed));
    }
    else
    {
        logger_->warn(
            "  ! REML did not converge ({} iterations in {:.2f}s)",
            max_iter,
            elapsed);
        logger_->warn(
            "    Try to increase max_iter or check the model specification.");
    }
    logger_->info("");
}

void RemlLogger::log_model_fit(const FreqModel& model, double loglike)
{
    logger_->info("  Model Fit:");
    logger_->info("  - AIC : {:.2f}", statistics::compute_aic(model, loglike));
    logger_->info("  - BIC : {:.2f}", statistics::compute_bic(model, loglike));
    logger_->info("");
}

void RemlLogger::log_fixed_effects(
    const FreqModel& model,
    const FreqState& state)
{
    logger_->info("  Fixed Effects:");
    logger_->info("  {:12} {:>12} {:>12}", "Effect", "Estimate", "SE");
    logger_->info(table_separator(40));
    for (Eigen::Index i = 0; i < state.fixed().coeff.size(); ++i)
    {
        std::string name = fmt::format("X{}", i);
        if (static_cast<size_t>(i) < model.fixed().names.size())
        {
            name = model.fixed().names[i];
        }
        logger_->info(
            "  {:12} {:>12.3f} {:>12.3f}",
            name,
            state.fixed().coeff(i),
            state.fixed().se(i));
    }
    logger_->info("");
}

void RemlLogger::log_variance_components(const FreqState& state)
{
    logger_->info("  Variance Components & Heritability:");
    logger_->info(
        "  {:12} {:>12} {:>12} {:>15} {:>12}",
        "Component",
        "Estimate",
        "SE",
        "Ratio (h²)",
        "SE");
    logger_->info(table_separator(69));

    // genetic effects with heritability
    double total_h2 = 0.0;
    for (const auto& g : state.genetic())
    {
        logger_->info(
            "  {:12} {:>12.3f} {:>12.3f} {:>15.3f} {:>12.3f}",
            g.type,
            g.variance,
            g.variance_se,
            g.heritability,
            g.heritability_se);
        total_h2 += g.heritability;
    }
    // random effects (if any)
    for (const auto& r : state.random())
    {
        logger_->info(
            "  {:12} {:>12.3f} {:>12.3f} {:>15} {:>12}",
            r.name,
            r.variance,
            r.variance_se,
            "-",
            "-");
    }
    // residual
    logger_->info(
        "  {:12} {:>12.3f} {:>12.3f} {:>15} {:>12}",
        "Residual",
        state.residual().variance,
        state.residual().variance_se,
        "-",
        "-");

    // total genetic variance and heritability
    if (state.genetic().size() > 1)
    {
        double total_vg = 0.0;
        for (const auto& g : state.genetic())
        {
            total_vg += g.variance;
        }
        logger_->info(table_separator(69));
        logger_->info(
            "  {:12} {:>12.3f} {:>12} {:>15.3f} {:>12}",
            "Total Vg",
            total_vg,
            "-",
            total_h2,
            "-");
    }
}

}  // namespace detail
}  // namespace gelex
