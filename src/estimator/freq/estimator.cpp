#include "gelex/estimator/freq/estimator.h"

#include <chrono>

#include <fmt/format.h>

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
    bool verbose) -> void
{
    auto start = std::chrono::steady_clock::now();

    if (!verbose)
    {
        logger_->set_level(spdlog::level::warn);
    }

    // log model info
    logger_->info(title(" GBLUP MODEL ANALYSIS "));
    logger_->info(subtitle("Model Specification"));
    logger_->info(item("Samples:  {:d}", model.num_individuals()));
    logger_->info(item("Fixed effects:  {:d}", model.fixed().X.cols()));
    logger_->info(item("Random effects:  {:d}", model.random().size()));
    logger_->info(item("Genetic effects:  {:d}", model.genetic().size()));
    logger_->info("");

    logger_->info(subtitle("Optimizer Specification"));
    logger_->info(item("Method:  {}", cyan("AI")));
    logger_->info(item("tolerance:  {:.2e}", tol_));
    logger_->info(item("Max Iterations:  {:d}", max_iter_));
    logger_->info("");

    OptimizerState opt_state(model);

    // EM initialization
    if (em_init)
    {
        em_step(model, state, opt_state);
    }

    // log iteration header
    logger_->info(title(" REML ESTIMATION "));
    std::string var_header;
    var_header += fmt::format("{:>12}", "V[e]");
    for (const auto& r : state.random())
    {
        var_header += fmt::format("{:>12}", fmt::format("V[{}]", r.name));
    }
    for (const auto& g : state.genetic())
    {
        var_header += fmt::format("{:>12}", fmt::format("V[{}]", g.name));
    }
    logger_->info(
        "{:>6} {:>12} {} {:>10}", "Iter.", "logL", var_header, "time");

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
        var_values += fmt::format("{:>12.4f}", state.residual().variance);
        for (const auto& r : state.random())
        {
            var_values += fmt::format("{:>12.4f}", r.variance);
        }
        for (const auto& g : state.genetic())
        {
            var_values += fmt::format("{:>12.4f}", g.variance);
        }
        logger_->info(
            "{:>6} {:>12.4f} {} {:>9.3f}s",
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
    logger_->info("Initializing with {} algorithm", cyan("EM"));

    std::string var_values;
    var_values += fmt::format("{:.4f}", state.residual().variance);
    for (const auto& r : state.random())
    {
        var_values += fmt::format(", {:.4f}", r.variance);
    }
    for (const auto& g : state.genetic())
    {
        var_values += fmt::format(", {:.4f}", g.variance);
    }

    logger_->info(
        "Initial: logL={:.4f} | \u03C3\u00B2=[{}] ({:.3f}s)",
        loglike,
        rebecca_purple(var_values),
        time_cost);
}

auto Estimator::report_results(
    const FreqModel& model,
    const FreqState& state,
    const OptimizerState& /*opt_state*/,
    double elapsed) -> void
{
    logger_->info(title(" RESULT "));

    // convergence status
    logger_->info(subtitle("Convergence"));
    if (converged_)
    {
        logger_->info(
            " \u25AA Status:  {} ({} iterations in {:.3f}s)",
            green("Success"),
            iter_count_,
            elapsed);
    }
    else
    {
        logger_->warn(
            " \u25AA Status:  {} ({} iterations in {:.3f}s)",
            red("Failed"),
            max_iter_,
            elapsed);
        logger_->warn(
            "Try to increase the max_iter or check the model specification.");
    }
    logger_->info(
        " \u25AA AIC:  {:.4f}", statistics::compute_aic(model, loglike_));
    logger_->info(
        " \u25AA BIC:  {:.4f}", statistics::compute_bic(model, loglike_));
    logger_->info("");

    // fixed effects
    logger_->info(subtitle("Fixed Effects"));
    for (Eigen::Index i = 0; i < state.fixed().coeff.size(); ++i)
    {
        std::string name = fmt::format("X{}", i);
        if (static_cast<size_t>(i) < model.fixed().names.size())
        {
            name = model.fixed().names[i];
        }
        logger_->info(item(
            "{}:  {}",
            name,
            with_std(state.fixed().coeff(i), state.fixed().se(i))));
    }
    logger_->info("");

    // variance components
    logger_->info(subtitle("Variance Components"));
    logger_->info(item(
        "Residual:  {}",
        with_std(state.residual().variance, state.residual().variance_se)));
    for (const auto& r : state.random())
    {
        logger_->info(
            item("{}:  {}", r.name, with_std(r.variance, r.variance_se)));
    }
    for (const auto& g : state.genetic())
    {
        logger_->info(
            item("{}:  {}", g.name, with_std(g.variance, g.variance_se)));
    }
    logger_->info("");

    // heritability
    if (!state.genetic().empty())
    {
        logger_->info(subtitle("Heritability"));
        for (const auto& g : state.genetic())
        {
            logger_->info(item(
                "{}:  {}",
                g.name,
                with_std(g.heritability, g.heritability_se)));
        }
    }

    logger_->info(title(""));
}

}  // namespace gelex
