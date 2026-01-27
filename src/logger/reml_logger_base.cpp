#include "reml_logger_base.h"

namespace gelex::detail
{

void RemlLoggerBase::set_verbose(bool verbose)
{
    if (!verbose)
    {
        logger_->set_level(spdlog::level::warn);
    }
}

void RemlLoggerBase::log_em_init(const FreqState& /*state*/, double /*loglike*/)
{
}

void RemlLoggerBase::log_iter_header(const FreqState& /*state*/) {}

void RemlLoggerBase::log_iteration(
    size_t /*iter*/,
    double /*loglike*/,
    const FreqState& /*state*/,
    double /*time_cost*/)
{
}

void RemlLoggerBase::log_iter_footer() {}

void RemlLoggerBase::log_results(
    const FreqModel& /*model*/,
    const FreqState& /*state*/,
    double /*loglike*/,
    bool /*converged*/,
    size_t /*iter_count*/,
    size_t /*max_iter*/,
    double /*elapsed*/)
{
}

}  // namespace gelex::detail
