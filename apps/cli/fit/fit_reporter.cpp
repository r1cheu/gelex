/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "fit_reporter.h"

#include <iterator>

#include <fmt/format.h>

#include "config.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/utils/formatter.h"
#include "gelex/model/bayes/model.h"
#include "gelex/types/mcmc_results.h"

namespace gelex::cli
{

namespace
{
const int kTableWidth = 40;
}  // namespace

FitReporter::FitReporter() : logger_(gelex::logging::get()) {}

auto FitReporter::on_event(const FitConfigLoadedEvent& event) const -> void
{
    logger_->info(
        gelex::command_banner(PROJECT_VERSION, "Model Fitting (MCMC)"));
    logger_->info("");
    logger_->info(gelex::section("[Config]"));
    logger_->info("  {:<12}: {}", "Method", fmt::format("{}", event.method));
    logger_->info(
        "  {:<12}: {} iters ({} burn-in, {} sampling)",
        "Chain",
        event.n_iters,
        event.n_burnin,
        event.n_iters - event.n_burnin);
    logger_->info("  {:<12}: {}", "Seed", event.seed);
    logger_->info("");
}

auto FitReporter::on_event(const FitModelReadyEvent& event) const -> void
{
    logger_->info(gelex::section("[Model Configuration]"));

    for (const auto& effect : event.model->random())
    {
        print_random_prior(effect);
    }

    print_genetic_prior(event.model->additive(), "Additive effect:");
    print_genetic_prior(event.model->dominant(), "Dominance effect:");
    print_residual_prior(event.model->residual());
}

auto FitReporter::on_event(const FitMcmcProgressEvent& event) -> void
{
    if (!init_progress_)
    {
        init_progress_ = true;
        logger_->info("");
        logger_->info(gelex::section("[MCMC Sampling]"));
        bar_ = detail::create_progress_bar(
            iter_, event.total, "{bar} {value}/{total} [{speed:.1f}/s]");
        bar_.display->show();
    }

    if (event.done)
    {
        bar_.display->done();
        logger_->info("");
        return;
    }

    iter_ = event.current;

    stats_.clear();
    if (event.h2)
    {
        fmt::format_to(std::back_inserter(stats_), "h²: {:.3f}", *event.h2);
    }
    if (event.d2)
    {
        fmt::format_to(
            std::back_inserter(stats_),
            "{}δ²: {:.3f}",
            stats_.empty() ? "" : " | ",
            *event.d2);
    }
    if (event.sigma2_e)
    {
        fmt::format_to(
            std::back_inserter(stats_),
            "{}σ²_e: {:.3f}",
            stats_.empty() ? "" : " | ",
            *event.sigma2_e);
    }

    if (bar_.after_bar)
    {
        bar_.after_bar->message(stats_);
    }
}

auto FitReporter::on_event(const FitMcmcCompleteEvent& event) const -> void
{
    print_fixed_summary(*event.result, *event.model, event.samples_collected);
    print_genetic_summary(
        event.result->additive(),
        event.model->additive(),
        GeneticEffectType::Add);
    print_genetic_summary(
        event.result->dominant(),
        event.model->dominant(),
        GeneticEffectType::Dom);
    print_residual_summary(*event.result);
}

auto FitReporter::on_event(const FitResultsSavedEvent& event) const -> void
{
    logger_->info(
        gelex::success(
            "Results saved to '{}' (.param, .snp.eff, .log)",
            event.out_prefix));
}

// --- Private helpers ---

auto FitReporter::print_variance_prior(
    const detail::ScaledInvChiSqParams& prior,
    double init_variance) const -> void
{
    logger_->info(
        gelex::subtask(
            "Variance: Scaled Inv-χ²(ν={:.4f}, S²={:.4f}), init: {:.4f}",
            prior.nu,
            prior.s2,
            init_variance));
}

auto FitReporter::print_summary_row(
    std::string_view name,
    const PosteriorSummary& summary,
    Eigen::Index index) const -> void
{
    logger_->info(
        "  {:<8} {:>10.6f} {:>10.6f}",
        name,
        summary.mean(index),
        summary.stddev(index));
}

auto FitReporter::print_random_prior(const bayes::RandomEffect& effect) const
    -> void
{
    std::string name = effect.levels ? effect.levels.value()[0] : "intercept";
    logger_->info(gelex::task("{}(rand)", name));
    print_variance_prior(effect.prior, effect.init_variance);
}

auto FitReporter::print_genetic_prior(
    const bayes::GeneticEffect* effect,
    std::string_view label) const -> void
{
    if (effect == nullptr)
    {
        return;
    }
    logger_->info(gelex::task("{}", label));
    print_variance_prior(
        effect->marker_variance_prior, effect->init_marker_variance);

    if (effect->init_pi && effect->init_pi->size() > 1)
    {
        std::string pi_str = "[";
        for (Eigen::Index i = 0; i < effect->init_pi->size(); ++i)
        {
            if (i > 0)
            {
                pi_str += ", ";
            }
            pi_str += fmt::format("{:.3f}", (*effect->init_pi)(i));
        }
        pi_str += "]";
        logger_->info(gelex::subtask("Mixture: {}", pi_str));
    }
}

auto FitReporter::print_residual_prior(const bayes::Residual& residual) const
    -> void
{
    logger_->info(gelex::task("Residual:"));
    print_variance_prior(residual.prior, residual.init_variance);
}

auto FitReporter::print_fixed_summary(
    const MCMCResult& result,
    const BayesModel& model,
    std::ptrdiff_t samples_collected) const -> void
{
    auto [effect, res] = std::make_pair(model.fixed(), result.fixed());
    if ((effect == nullptr) || (res == nullptr))
    {
        return;
    }

    logger_->info("");
    logger_->info(gelex::section("[Posterior Summary]"));
    logger_->info("  Samples collected per parameter: {}", samples_collected);
    logger_->info("");

    logger_->info("  {:<8} {:>8} {:>8}", "Parameter", "Mean", "SD");
    logger_->info(gelex::table_separator(kTableWidth));

    Eigen::Index idx{};
    for (const auto& covariates : effect->levels)
    {
        if (covariates)
        {
            for (const auto& level : covariates.value())
            {
                print_summary_row(level, res->coeffs, idx);
                ++idx;
            }
        }
        else
        {
            print_summary_row(effect->names[idx], res->coeffs, idx);
            ++idx;
        }
    }
}

auto FitReporter::print_genetic_summary(
    const BaseMarkerSummary* summary,
    const bayes::GeneticEffect* effect,
    GeneticEffectType type) const -> void
{
    if (summary == nullptr)
    {
        return;
    }

    std::string label
        = type == GeneticEffectType::Dom ? "Dominance" : "Additive";
    std::string h_name = type == GeneticEffectType::Dom ? "δ²" : "h²";

    logger_->info(gelex::named_section(label, kTableWidth, 2));
    print_summary_row("σ²", summary->variance);
    print_summary_row(h_name, summary->heritability);

    if ((effect != nullptr) && effect->init_pi && effect->init_pi->size() > 1)
    {
        for (Eigen::Index i = 0; i < summary->mixture_proportion.size(); ++i)
        {
            print_summary_row(
                fmt::format("π[{}]", i), summary->mixture_proportion, i);
        }
        for (Eigen::Index i = 0; i < summary->component_variance.size(); ++i)
        {
            print_summary_row(
                fmt::format("σ²[{}]", i + 1), summary->component_variance, i);
        }
    }
}

auto FitReporter::print_residual_summary(const MCMCResult& result) const -> void
{
    logger_->info(gelex::named_section("Residual", kTableWidth, 2));
    print_summary_row("σ²", result.residual());
    logger_->info(gelex::table_separator(kTableWidth));
    logger_->info("");
}

}  // namespace gelex::cli
