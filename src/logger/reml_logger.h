#ifndef GELEX_LOGGER_REML_LOGGER_H_
#define GELEX_LOGGER_REML_LOGGER_H_

#include "reml_logger_base.h"

#include <cstddef>

namespace gelex
{
class FreqModel;
class FreqState;

namespace detail
{

class RemlLogger : public RemlLoggerBase
{
   public:
    RemlLogger() = default;

    void set_verbose(bool verbose) override;
    void log_em_init(const FreqState& state, double loglike) override;
    void log_iter_header(const FreqState& state) override;
    void log_iteration(
        size_t iter,
        double loglike,
        const FreqState& state,
        double time_cost) override;
    void log_iter_footer() override;
    void log_results(
        const FreqModel& model,
        const FreqState& state,
        double loglike,
        bool converged,
        size_t iter_count,
        size_t max_iter,
        double elapsed) override;

   private:
    void log_convergence(
        bool converged,
        size_t iter_count,
        size_t max_iter,
        double elapsed);
    void log_model_fit(const FreqModel& model, double loglike);
    void log_fixed_effects(const FreqModel& model, const FreqState& state);
    void log_variance_components(const FreqState& state);
};

}  // namespace detail
}  // namespace gelex

#endif  // GELEX_LOGGER_REML_LOGGER_H_
