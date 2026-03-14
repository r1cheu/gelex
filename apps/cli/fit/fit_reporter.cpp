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
#include "gelex/infra/logging/formatter.h"
#include "gelex/model/bayes/model.h"
#include "gelex/types/genetic_effect_type.h"
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

    for (const auto& effect : event.model->genetics())
    {
        print_genetic_prior(
            &effect,
            fmt::format("{} effect:", genetic_kind::to_string(effect.type)));
    }
    print_residual_prior(event.model->residual());
}

auto FitReporter::on_event(const FitMcmcProgressEvent& event) -> void
{
    if (!init_progress_)
    {
        init_progress_ = true;
        logger_->info("");
        logger_->info(gelex::section("[MCMC Sampling]"));
        bar_ = create_progress_bar(
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
    for (const auto& [type, h] : event.genetic_heritabilities)
    {
        fmt::format_to(
            std::back_inserter(stats_),
            "{}{}:{:.3f}",
            stats_.empty() ? "" : " | ",
            genetic_kind::to_heritability_label(type),
            h);
    }
    if (event.sigma2_e)
    {
        fmt::format_to(
            std::back_inserter(stats_),
            "{}σ²_e: {:.3f}",
            stats_.empty() ? "" : " | ",
            *event.sigma2_e);
    }
    if (event.dom_positive_prob)
    {
        fmt::format_to(
            std::back_inserter(stats_),
            " | p: {:.3f}",
            *event.dom_positive_prob);
    }

    if (bar_.after_bar)
    {
        bar_.after_bar->message(stats_);
    }
}

auto FitReporter::on_event(const FitMcmcCompleteEvent& event) const -> void
{
    print_fixed_summary(*event.result, event.samples_collected);
    for (const auto& summary : event.result->genetics())
    {
        print_genetic_summary(
            &summary, event.model->genetic(summary.type), summary.type);
    }
    print_residual_summary(*event.result);
}

auto FitReporter::on_event(const FitResultsSavedEvent& event) const -> void
{
    logger_->info(
        gelex::success(
            "Results saved to '{}' (.param, .summary, .snp.eff, .log)",
            event.out_prefix));
}

// --- Private helpers ---

auto FitReporter::print_variance_prior(
    const detail::ScaledInvChiSqParams& prior,
    double init_variance) const -> void
{
    logger_->info(
        "    Variance: Scaled Inv-χ²(ν={:.4f}, S²={:.4f}), init: {:.4f}",
        prior.nu,
        prior.s2,
        init_variance);
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
    logger_->info("   {}(rand)", name);
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
    logger_->info("   {}", label);
    print_variance_prior(
        effect->marker_variance_prior, effect->init_marker_variance);

    if (effect->mixture && effect->mixture->init_proportion.size() > 1)
    {
        const auto& init_proportion = effect->mixture->init_proportion;
        std::string pi_str = "[";
        for (Eigen::Index i = 0; i < init_proportion.size(); ++i)
        {
            if (i > 0)
            {
                pi_str += ", ";
            }
            pi_str += fmt::format("{:.3f}", init_proportion(i));
        }
        pi_str += "]";
        logger_->info("    Mixture: {}", pi_str);
    }
}

auto FitReporter::print_residual_prior(const bayes::Residual& residual) const
    -> void
{
    logger_->info("   Residual:");
    print_variance_prior(residual.prior, residual.init_variance);
}

auto FitReporter::print_fixed_summary(
    const MCMCResult& result,
    std::ptrdiff_t samples_collected) const -> void
{
    const auto& fixed = result.fixed();

    logger_->info("");
    logger_->info(gelex::section("[Posterior Summary]"));
    logger_->info("  Samples collected per parameter: {}", samples_collected);
    logger_->info("");

    logger_->info("  {:<8} {:>8} {:>8}", "Parameter", "Mean", "SD");
    logger_->info(gelex::table_separator(kTableWidth));

    fixed.for_each_term([&](const std::string& term, Eigen::Index i)
                        { print_summary_row(term, fixed.coeffs, i); });
}

auto FitReporter::print_genetic_summary(
    const GeneticSummary* summary,
    const bayes::GeneticEffect* effect,
    GeneticKind type) const -> void
{
    if (summary == nullptr)
    {
        return;
    }

    std::string label{genetic_kind::to_string(type)};
    std::string h_name{genetic_kind::to_heritability_label(type)};

    logger_->info(gelex::named_section(label, kTableWidth, 2));
    print_summary_row("σ²", summary->variance);
    print_summary_row(h_name, summary->heritability);

    if (summary->mixture)
    {
        const auto& mix = *summary->mixture;
        for (Eigen::Index i = 0; i < mix.mixture_proportion.size(); ++i)
        {
            print_summary_row(
                fmt::format("π[{}]", i), mix.mixture_proportion, i);
        }
        for (Eigen::Index i = 0; i < mix.component_variance.size(); ++i)
        {
            print_summary_row(
                fmt::format("σ²[{}]", i + 1), mix.component_variance, i);
        }
    }

    if (summary->sign)
    {
        print_summary_row("p(+)", summary->sign->positive_prob);
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
