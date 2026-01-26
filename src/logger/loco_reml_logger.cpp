#include "loco_reml_logger.h"

#include <fmt/format.h>

#include "../utils/formatter.h"
#include "gelex/model/freq/model.h"

namespace gelex::detail
{

LocoRemlLogger::LocoRemlLogger(std::string chr_name)
    : chr_name_(std::move(chr_name))
{
}

void LocoRemlLogger::set_verbose(bool /*verbose*/)
{
    // LOCO mode keeps info level to show concise output
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
    std::string h2_str;
    for (const auto& g : state.genetic())
    {
        if (!h2_str.empty())
        {
            h2_str += ", ";
        }
        h2_str += fmt::format("h²({})={:.2f}", g.type, g.heritability);
    }

    std::string status = converged ? "" : " [!]";
    logger_->info(
        "   {} Chr {:>2}: LogL={:>8.2f} | {} | Time={:.2f}s{}",
        fmt::format(fmt::fg(fmt::color::light_cyan), "●"),
        chr_name_,
        loglike,
        h2_str,
        elapsed,
        status);
}

}  // namespace gelex::detail
