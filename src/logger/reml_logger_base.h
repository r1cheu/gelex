#ifndef GELEX_LOGGER_REML_LOGGER_BASE_H_
#define GELEX_LOGGER_REML_LOGGER_BASE_H_

#include <cstddef>
#include <memory>

#include <spdlog/logger.h>

namespace gelex
{
class FreqModel;
class FreqState;

namespace detail
{

class RemlLoggerBase
{
   public:
    RemlLoggerBase();
    virtual ~RemlLoggerBase() = default;

    virtual void set_verbose(bool verbose);
    virtual void log_em_init(const FreqState& state, double loglike);
    virtual void log_iter_header(const FreqState& state);
    virtual void log_iteration(
        size_t iter,
        double loglike,
        const FreqState& state,
        double time_cost);
    virtual void log_iter_footer();
    virtual void log_results(
        const FreqModel& model,
        const FreqState& state,
        double loglike,
        bool converged,
        size_t iter_count,
        size_t max_iter,
        double elapsed);

   protected:
    std::shared_ptr<spdlog::logger> logger_;
};

}  // namespace detail
}  // namespace gelex

#endif  // GELEX_LOGGER_REML_LOGGER_BASE_H_
