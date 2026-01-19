#include "gelex/estimator/freq/estimator.h"

#include <Eigen/Core>
#include <chrono>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include "../src/utils/formatter.h"
#include "../src/utils/utils.h"
#include "gelex/estimator/freq/effect_solver.h"
#include "gelex/estimator/freq/statistics.h"
#include "gelex/logger.h"
#include "gelex/model/freq/model.h"
#include "gelex/optim/policy.h"
#include "gelex/optim/variance_calculator.h"

namespace gelex
{

Estimator::Estimator(size_t max_iter, double tol)
    : optimizer_(tol), max_iter_(max_iter), tol_(tol), logger_(logging::get())
{
}

auto Estimator::fit(
    const FreqModel& model,
    FreqState& state,
    bool em_init,
    bool verbose) -> Eigen::MatrixXd
{
    auto start = std::chrono::steady_clock::now();

    if (!verbose)
    {
        logger_->set_level(spdlog::level::warn);
    }

    OptimizerState opt_state(model);

    // EM initialization
    if (em_init)
    {
        em_step(model, state, opt_state);
    }

    // log iteration header
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

    // AI iterations
    for (size_t iter = 1; iter <= max_iter_; ++iter)
    {
        double time_cost{};
        {
            Timer timer{time_cost};
            optimizer_.step<AIPolicy>(model, state, opt_state);
        }

        loglike_ = variance_calculator::compute_loglike(model, opt_state);

        // log iteration
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
            "  {:<4} {:>12.2f} {} {:>9.2f}s",
            iter,
            loglike_,
            var_values,
            time_cost);

        if (optimizer_.is_converged())
        {
            converged_ = true;
            iter_count_ = iter;
            break;
        }
    }
    logger_->info(gelex::table_separator(67));

    if (!converged_)
    {
        iter_count_ = max_iter_;
    }

    // compute final results
    effect_solver::compute_fixed_effects(model, state, opt_state);
    effect_solver::compute_random_effects(model, state, opt_state);
    statistics::compute_variance_se(state, opt_state);
    statistics::compute_heritability(state, opt_state);

    auto elapsed = std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - start)
                       .count();
    report_results(model, state, opt_state, elapsed);
    return std::move(opt_state.v);
}

auto Estimator::em_step(
    const FreqModel& model,
    FreqState& state,
    OptimizerState& opt_state) -> void
{
    double time_cost{};
    {
        Timer timer{time_cost};
        optimizer_.step<EMPolicy>(model, state, opt_state);
    }

    double loglike = variance_calculator::compute_loglike(model, opt_state);

    // log EM initialization
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

auto Estimator::report_results(
    const FreqModel& model,
    const FreqState& state,
    const OptimizerState& /*opt_state*/,
    double elapsed) -> void
{
    logger_->info("");
    logger_->info(
        fmt::format(
            fmt::fg(fmt::color::light_cyan),
            "── REML Results {}",
            separator(70 - 16)));

    // convergence status
    if (converged_)
    {
        logger_->info(success(
            "Converged successfully in {} iterations ({:.2f}s)",
            iter_count_,
            elapsed));
    }
    else
    {
        logger_->warn(
            "  ! REML did not converge ({} iterations in {:.2f}s)",
            max_iter_,
            elapsed);
        logger_->warn(
            "    Try to increase max_iter or check the model specification.");
    }
    logger_->info("");

    logger_->info("  Model Fit:");
    logger_->info("  - AIC : {:.2f}", statistics::compute_aic(model, loglike_));
    logger_->info("  - BIC : {:.2f}", statistics::compute_bic(model, loglike_));
    logger_->info("");

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

    // Variance components & heritability table
    logger_->info("  Variance Components & Heritability:");
    logger_->info(
        "  {:12} {:>12} {:>12} {:>15} {:>12}",
        "Component",
        "Estimate",
        "SE",
        "Ratio (h²)",
        "SE");
    logger_->info(table_separator(69));

    // Genetic effects with heritability
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
    // Random effects (if any)
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
    // Residual
    logger_->info(
        "  {:12} {:>12.3f} {:>12.3f} {:>15} {:>12}",
        "Residual",
        state.residual().variance,
        state.residual().variance_se,
        "-",
        "-");

    // Total genetic variance and heritability
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

    logger_->info(
        fmt::format(fmt::fg(fmt::color::light_cyan), "{}", separator(70)));
}

}  // namespace gelex
