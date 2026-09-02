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

#include "reporter.h"

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <fmt/base.h>
#include <fmt/format.h>
#include <iterator>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <variant>

#include "gelex/algo/mcmc/result.h"
#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/legacy_genetic_prior.h"
#include "gelex/bayes/labels.h"
#include "gelex/bayes/legacy_prior.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/parameter/values.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/stats/result.h"
#include "gelex/types/genetic_mode.h"

#include "cli/formatter.h"
#include "cli/progress_bar.h"
#include "cli/report_printer.h"
#include "cli/table.h"

namespace cli
{

GenoReporter::GenoReporter() : eta_(1) {}

auto GenoReporter::show_total(int64_t num_variants) const -> void
{
    cli::printer().line(
        "{}", cli::field("Variants", "{}", cli::AbbrNumber(num_variants)));
}

auto GenoReporter::show_loaded(const gelex::bayes::GeneticDesign& design) const
    -> void
{
    for (const gelex::GeneticMode mode : design.modes().each())
    {
        show_loaded_mode(
            mode,
            design.cols(),
            design.cols()
                - static_cast<Eigen::Index>(design.valid_indices(mode).size()));
    }
}

auto GenoReporter::show_loaded_mode(
    gelex::GeneticMode mode,
    int64_t num_snps,
    int64_t invalid_snps) const -> void
{
    const std::string label
        = (mode == gelex::GeneticMode::D) ? "Dominance" : "Additive";
    if (invalid_snps == 0)
    {
        cli::printer().line("{}", cli::field(label, "all valid"));
        return;
    }
    cli::printer().line(
        "{}",
        cli::field(
            label,
            "{} valid ({} excluded)",
            cli::AbbrNumber(num_snps - invalid_snps),
            cli::AbbrNumber(invalid_snps)));
}

auto GenoReporter::on_event(const gelex::GenotypeProgressEvent& event) -> void
{
    if (event.done)
    {
        if (bar_active_)
        {
            bar_.display->done();
            bar_active_ = false;
            cli::clear_finished_line();
            cli::printer().on_progress_finished();
        }
        return;
    }

    if (!bar_active_)
    {
        total_ = event.total;
        progress_ = 0;
        eta_.reset(total_);
        bar_ = cli::create_progress_bar(progress_, total_);
        bar_.display->show();
        bar_active_ = true;
    }

    progress_ = event.current;
    if (bar_.after_bar)
    {
        bar_.after_bar->message(
            fmt::format(
                "{:.1f}% ({}/{} SNPs) | ETA: {}",
                static_cast<double>(progress_) / static_cast<double>(total_)
                    * 100.0,
                cli::AbbrNumber(progress_),
                cli::AbbrNumber(total_),
                eta_.get_eta(progress_)));
    }
}

auto McmcReporter::show_dataset_summary(
    const gelex::BayesModel& model,
    std::string_view pheno_name) -> void
{
    auto& p = cli::printer();
    p.block(cli::section("Dataset Summary:"));
    p.line(cli::field("Trait", "{}", pheno_name));
    p.line(cli::field("Analyzed Samples", "{}", model.num_individuals()));
    p.line(cli::field("Covariates", "{}", model.fixed().X.cols()));
}

auto McmcReporter::show_prior(
    const gelex::bayes::BayesPrior& prior,
    const gelex::BayesModel& model) -> void
{
    prior_ = &prior;

    cli::printer().block(cli::section("Prior Configuration:"));
    if (!model.random().empty())
    {
        print_random_prior(prior.random());
    }
    for (const auto& genetic : prior.genetics())
    {
        print_genetic_prior(genetic);
    }
    print_residual_prior(prior.residual());
}

auto McmcReporter::print_variance_prior(
    const gelex::bayes::ScaledInvChiSqPrior& prior,
    double init_variance) -> void
{
    cli::printer().line(
        "    Variance: Scaled Inv-χ²(ν={:.4f}, S²={:.4g}), init: {:.4g}",
        prior.degrees_of_freedom(),
        prior.scale(),
        init_variance);
}

auto McmcReporter::print_random_prior(const gelex::bayes::RandomPrior& prior)
    -> void
{
    cli::printer().line("   Random effect:");
    print_variance_prior(prior.prior(), prior.initial_value());
}

auto McmcReporter::print_genetic_prior(const gelex::bayes::GeneticPrior& prior)
    -> void
{
    auto format_vec = [](const auto& p)
    {
        auto formatted = std::span(p.data(), p.size())
                         | std::views::transform(
                             [](double v) { return fmt::format("{:.3f}", v); });
        return fmt::join(formatted, ", ");
    };

    auto print_proportion = [&format_vec](const auto& proportion)
    {
        cli::printer().line(
            "    Proportion: [{}]", format_vec(proportion.initial_value()));
        cli::printer().line(
            "    Proportion update: {}",
            proportion.is_sampled() ? "yes" : "no");
    };

    std::visit(
        [&](const auto& genetic_group)
        {
            using Group = std::remove_cvref_t<decltype(genetic_group)>;
            if constexpr (
                std::is_same_v<Group, gelex::bayes::SingleGeneticPrior>)
            {
                std::visit(
                    [&](const auto& genetic)
                    {
                        cli::printer().line("   {} effect:", genetic.mode());
                        const auto& parameter = genetic.variance().parameter();
                        print_variance_prior(
                            parameter.prior(), parameter.initial_value());
                        if constexpr (requires { genetic.proportion(); })
                        {
                            print_proportion(genetic.proportion());
                        }
                        if constexpr (requires { genetic.multiplier(); })
                        {
                            cli::printer().line(
                                "    Multiplier: [{}]",
                                format_vec(genetic.multiplier()));
                        }
                    },
                    genetic_group);
            }
            else
            {
                cli::printer().line("   A, D effect:");
                std::visit(
                    [&](const auto& genetic)
                    {
                        for (const auto mode :
                             {gelex::GeneticMode::A, gelex::GeneticMode::D})
                        {
                            const auto& parameter
                                = genetic.variance(mode).parameter();
                            print_variance_prior(
                                parameter.prior(), parameter.initial_value());
                        }
                        print_proportion(genetic.proportion());
                    },
                    genetic_group);
            }
        },
        prior);
}

auto McmcReporter::print_residual_prior(
    const gelex::bayes::ResidualPrior& prior) -> void
{
    cli::printer().line("   Residual:");
    print_variance_prior(prior.prior(), prior.initial_value());
}

auto McmcReporter::on_event(const gelex::MCMCProgressEvent& event) -> void
{
    if (!init_progress_)
    {
        init_progress_ = true;
        cli::printer().block(cli::section("MCMC Sampling:"));
        bar_ = cli::create_progress_bar(
            iter_, event.total, "{bar} {value}/{total} [{speed:.1f}/s]");
        bar_.display->show();
    }

    if (event.done)
    {
        bar_.display->done();
        cli::printer().on_progress_finished();
        return;
    }

    iter_ = event.current;
    bar_.before_bar->message("");

    const auto* state = event.state;
    stats_.clear();
    if (prior_ != nullptr)
    {
        const auto priors = prior_->genetics();
        for (auto [i, gen] : std::views::enumerate(state->genetics()))
        {
            const auto& prior = priors[static_cast<std::size_t>(i)];
            std::visit(
                [&](const auto& prior_block, const auto& block)
                {
                    using Prior = std::decay_t<decltype(prior_block)>;
                    using Block = std::decay_t<decltype(block)>;
                    if constexpr (
                        std::is_same_v<Prior, gelex::bayes::SingleGeneticPrior>
                        && std::is_same_v<
                            Block,
                            gelex::bayes::SingleGeneticBlockState>)
                    {
                        const auto mode = gelex::bayes::mode(prior_block);
                        const auto& genetic = block.state();
                        fmt::format_to(
                            std::back_inserter(stats_),
                            "{}{}:{:.3f}",
                            stats_.empty() ? "" : " | ",
                            gelex::bayes::to_heritability_label(mode),
                            genetic.heritability);
                        std::visit(
                            [&](const auto& prior_state)
                            {
                                using PriorState
                                    = std::decay_t<decltype(prior_state)>;
                                if constexpr (
                                    std::is_same_v<
                                        PriorState,
                                        gelex::bayes::
                                            SingleSharedSpikeSlabGaussianState>
                                    || std::is_same_v<
                                        PriorState,
                                        gelex::bayes::
                                            SinglePerMarkerSpikeSlabGaussianState>
                                    || std::is_same_v<
                                        PriorState,
                                        gelex::bayes::
                                            SingleScaledMixtureGaussianState>)
                                {
                                    fmt::format_to(
                                        std::back_inserter(stats_), " | π:[");
                                    const auto& proportion
                                        = prior_state.proportion();
                                    for (Eigen::Index j = 0;
                                         j < proportion.size();
                                         ++j)
                                    {
                                        fmt::format_to(
                                            std::back_inserter(stats_),
                                            "{}{:.3f}",
                                            j == 0 ? "" : ",",
                                            proportion(j));
                                    }
                                    fmt::format_to(
                                        std::back_inserter(stats_), "]");
                                }
                            },
                            block.prior_state());
                    }
                    else if constexpr (
                        std::is_same_v<Prior, gelex::bayes::JointGeneticPrior>
                        && std::is_same_v<
                            Block,
                            gelex::bayes::JointGeneticBlockState>)
                    {
                        for (const auto mode :
                             {gelex::GeneticMode::A, gelex::GeneticMode::D})
                        {
                            const auto& genetic = block.state(mode);
                            fmt::format_to(
                                std::back_inserter(stats_),
                                "{}{}:{:.3f}",
                                stats_.empty() ? "" : " | ",
                                gelex::bayes::to_heritability_label(mode),
                                genetic.heritability);
                        }
                        std::visit(
                            [&](const auto& prior_state)
                            {
                                using PriorState
                                    = std::decay_t<decltype(prior_state)>;
                                if constexpr (
                                    std::is_same_v<
                                        PriorState,
                                        gelex::bayes::JointGaussianMixtureState>
                                    || std::is_same_v<
                                        PriorState,
                                        gelex::bayes::
                                            JointHalfNormalMixtureState>)
                                {
                                    fmt::format_to(
                                        std::back_inserter(stats_), " | π:[");
                                    const auto& proportion
                                        = prior_state.proportion();
                                    for (Eigen::Index j = 0;
                                         j < proportion.size();
                                         ++j)
                                    {
                                        fmt::format_to(
                                            std::back_inserter(stats_),
                                            "{}{:.3f}",
                                            j == 0 ? "" : ",",
                                            proportion(j));
                                    }
                                    fmt::format_to(
                                        std::back_inserter(stats_), "]");
                                }
                                if constexpr (
                                    std::is_same_v<
                                        PriorState,
                                        gelex::bayes::
                                            JointHalfNormalMixtureState>)
                                {
                                    fmt::format_to(
                                        std::back_inserter(stats_),
                                        " | p_s:{:.3f}",
                                        prior_state.dominance_sign()
                                            .positive_probability);
                                }
                            },
                            block.prior_state());
                    }
                },
                prior,
                gen);
        }
    }
    fmt::format_to(
        std::back_inserter(stats_),
        "{}σ²_e: {:.3f}",
        stats_.empty() ? "" : " | ",
        state->residual().variance);

    if (bar_.after_bar)
    {
        bar_.after_bar->message(stats_);
    }
}

namespace
{
auto effect_of(std::string_view path) -> std::string_view
{
    std::size_t start = 0;
    while (start < path.size())
    {
        const auto end = path.find('/', start);
        const auto segment = end == std::string_view::npos
                                 ? path.substr(start)
                                 : path.substr(start, end - start);
        if (segment == "A" || segment == "D")
        {
            return segment;
        }
        if (end == std::string_view::npos)
        {
            break;
        }
        start = end + 1;
    }
    return "-";
}
}  // namespace

auto McmcReporter::show_summary(const gelex::Result& result) -> void
{
    auto& p = cli::printer();
    p.block(cli::section("MCMC Summary:"));
    p.line("   Draws: {}", result.samples_collected());

    // Variance-scale parameters (σ², σ²_marker, σ²_e) span orders of magnitude,
    // so scientific notation avoids truncating tiny per-marker variances; the
    // bounded ratios (h²/δ², π, p_s) read best as fixed-point. Splitting keeps
    // each column single-scale and aligned. Classification is by report label:
    // variance labels alone start with "σ²".
    auto make_table = []
    {
        Table t;
        t.column("Parameter", Align::left);
        t.column("Effect", Align::right);
        t.column("Estimate", Align::right);
        t.column("SE", Align::right);
        return t;
    };
    Table variances = make_table();
    Table ratios = make_table();
    bool has_variance = false;
    bool has_ratio = false;

    for (const auto& record : result.records())
    {
        if (!std::holds_alternative<gelex::RunningStatsResult>(record.value)
            || !record.names
            || std::string_view{record.path}.ends_with("/coeffs"))
        {
            continue;
        }
        const auto& stats = std::get<gelex::RunningStatsResult>(record.value);
        const auto effect = std::string{effect_of(record.path)};
        for (const auto [i, name] : std::views::enumerate(*record.names))
        {
            if (std::string_view{name}.starts_with("σ²"))
            {
                variances.row(
                    {name,
                     effect,
                     fmt::format("{:.3e}", stats.mean(i)),
                     fmt::format("{:.3e}", stats.stddev(i))});
                has_variance = true;
            }
            else
            {
                ratios.row(
                    {name,
                     effect,
                     fmt::format("{:.4f}", stats.mean(i)),
                     fmt::format("{:.4f}", stats.stddev(i))});
                has_ratio = true;
            }
        }
    }

    if (has_variance)
    {
        p.ensure_blank();
        p.line("   Variance components:");
        p.line(variances.render());
    }
    if (has_ratio)
    {
        p.ensure_blank();
        p.line("   Ratios & probabilities:");
        p.line(ratios.render());
    }
}

auto McmcReporter::on_event(const gelex::MCMCCheckpointSavedEvent& /*event*/)
    -> void
{
    bar_.before_bar->message(fmt::format(" ckpt saved"));
}

}  // namespace cli
