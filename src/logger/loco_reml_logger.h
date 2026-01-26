#ifndef GELEX_LOGGER_LOCO_REML_LOGGER_H_
#define GELEX_LOGGER_LOCO_REML_LOGGER_H_

#include "reml_logger_base.h"

#include <string>

namespace gelex::detail
{

class LocoRemlLogger : public RemlLoggerBase
{
   public:
    explicit LocoRemlLogger(std::string chr_name);

    void set_verbose(bool verbose) override;
    void log_results(
        const FreqModel& model,
        const FreqState& state,
        double loglike,
        bool converged,
        size_t iter_count,
        size_t max_iter,
        double elapsed) override;

   private:
    std::string chr_name_;
};

}  // namespace gelex::detail

#endif  // GELEX_LOGGER_LOCO_REML_LOGGER_H_
