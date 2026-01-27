#include "loco_reml_logger.h"

#include <fmt/format.h>

#include "../utils/formatter.h"
#include "gelex/model/freq/model.h"

namespace gelex::detail
{

LocoRemlLogger::LocoRemlLogger(std::string chr_name)
    : chr_name_(std::move(chr_name))
{
    result_.chr_name = chr_name_;
}

void LocoRemlLogger::set_status(std::shared_ptr<barkeep::StatusDisplay> status)
{
    status_ = std::move(status);
}

void LocoRemlLogger::set_verbose(bool /*verbose*/)
{
    // LOCO mode keeps info level to show concise output
}

void LocoRemlLogger::log_iteration(
    size_t /*iter*/,
    double loglike,
    const FreqState& state,
    double time_cost)
{
    accumulated_time_ += time_cost;

    if (status_)
    {
        std::string h2_str;
        for (const auto& g : state.genetic())
        {
            if (!h2_str.empty())
            {
                h2_str += ", ";
            }
            h2_str += fmt::format("h²({})={:.2f}", g.type, g.heritability);
        }

        status_->message(
            fmt::format(
                "Processing Chr{} LogL={:.2f} {} Time={:.2f}s",
                chr_name_,
                loglike,
                h2_str,
                accumulated_time_));
    }
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
    std::string h2_str;
    for (const auto& g : state.genetic())
    {
        result_.genetic.push_back({
            .type = g.type,
            .variance = g.variance,
            .heritability = g.heritability,
        });

        if (!h2_str.empty())
        {
            h2_str += ", ";
        }
        h2_str += fmt::format("h²({})={:.2f}", g.type, g.heritability);
    }

    if (status_)
    {
        std::string status_text = converged ? "" : " [!]";
        status_->message(
            fmt::format(
                "Chr {:>2}: LogL={:>8.2f} | {} | Time={:.2f}s{}",
                chr_name_,
                loglike,
                h2_str,
                elapsed,
                status_text));
    }
}

}  // namespace gelex::detail
