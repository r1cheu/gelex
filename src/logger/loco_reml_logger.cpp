#include "loco_reml_logger.h"

#include <fmt/format.h>

#include "../utils/formatter.h"
#include "gelex/logger.h"
#include "gelex/model/freq/model.h"

namespace gelex::detail
{

LocoRemlLogger::LocoRemlLogger(std::string chr_name)
    : chr_name_(std::move(chr_name))
{
    result_.chr_name = chr_name_;
}

void LocoRemlLogger::set_verbose(bool /*verbose*/)
{
    // LOCO mode keeps info level to show concise output
}

void LocoRemlLogger::log_iteration(
    size_t /*iter*/,
    double /*loglike*/,
    const FreqState& /*state*/,
    double /*time_cost*/)
{
}

void LocoRemlLogger::log_results(
    const FreqModel& /*model*/,
    const FreqState& state,
    double loglike,
    bool converged,
    size_t /*iter_count*/,
    size_t /*max_iter*/,
    double elapsed)
{
    result_.loglike = loglike;
    result_.converged = converged;
    result_.elapsed = elapsed;
    result_.residual_variance = state.residual().variance;

    result_.genetic.clear();
    for (const auto& g : state.genetic())
    {
        result_.genetic.push_back({
            .type = g.type,
            .variance = g.variance,
            .heritability = g.heritability,
        });
    }
}

void print_loco_reml_summary(const std::vector<LocoRemlResult>& results)
{
    if (results.empty())
    {
        return;
    }

    auto logger = gelex::logging::get();
    logger->info("");

    size_t num_grm = results[0].genetic.size();
    // helper to format variance columns
    auto format_variances = [num_grm](const auto& values) -> std::string
    {
        std::string row;
        for (size_t i = 0; i < num_grm; ++i)
        {
            row += fmt::format("  {:>10.4f}", values[i]);
        }
        return row;
    };

    // build header
    std::string header = fmt::format("  {:>5}  {:>10}", "Chr", "LogL");
    for (const auto& g : results[0].genetic)
    {
        header += fmt::format("  {:>10}", fmt::format("V({})", g.type));
    }
    header += fmt::format("  {:>10}  {:>7}  {:>4}", "V(e)", "Time", "Conv");

    logger->info(
        fmt::format(
            fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
            "── LOCO REML Summary {}",
            separator(70 - 21)));
    logger->info("{}", header);
    logger->info("{}", table_separator());

    // accumulate statistics while printing rows
    std::vector<double> sum_vg(num_grm, 0.0);
    std::vector<double> sum_h2(num_grm, 0.0);
    double sum_ve = 0.0;

    for (const auto& r : results)
    {
        std::string conv_mark
            = r.converged ? fmt::format(fmt::fg(fmt::color::light_green), "✓")
                          : fmt::format(fmt::fg(fmt::color::orange_red), "✗");

        std::vector<double> variances(num_grm);
        for (size_t i = 0; i < num_grm; ++i)
        {
            variances[i] = (i < r.genetic.size()) ? r.genetic[i].variance : 0.0;
            sum_vg[i] += variances[i];
            sum_h2[i]
                += (i < r.genetic.size()) ? r.genetic[i].heritability : 0.0;
        }

        logger->info(
            "  {:>5}  {:>10.2f}{}  {:>10.4f}  {:>6.2f}s    {}",
            r.chr_name,
            r.loglike,
            format_variances(variances),
            r.residual_variance,
            r.elapsed,
            conv_mark);

        sum_ve += r.residual_variance;
    }

    logger->info("{}", table_separator());

    // print summary rows
    auto n = static_cast<double>(results.size());
    std::vector<double> mean_vg(num_grm);
    std::vector<double> mean_h2(num_grm);
    for (size_t i = 0; i < num_grm; ++i)
    {
        mean_vg[i] = sum_vg[i] / n;
        mean_h2[i] = sum_h2[i] / n;
    }

    logger->info(
        "  {:>5}  {:>10}{}  {:>10.4f}",
        "Mean",
        "",
        format_variances(mean_vg),
        sum_ve / n);
    logger->info("  {:>5}  {:>10}{}", "h²", "", format_variances(mean_h2));

    logger->info(
        fmt::format(
            fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
            "{}",
            separator()));
}

}  // namespace gelex::detail
