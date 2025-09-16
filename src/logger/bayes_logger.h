#pragma once
#include <memory>

#include <gelex/barkeep.h>
#include <spdlog/logger.h>

#include "gelex/estimator/bayes/params.h"
#include "gelex/estimator/bayes/result.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{

namespace detail
{
class MCMCLogger
{
   public:
    MCMCLogger();
    void set_verbose(bool verbose);

    template <typename... Args>
    void info(spdlog::format_string_t<Args...> fmt, Args&&... args)
    {
        logger_->info(fmt, std::forward<Args>(args)...);
    }
    template <typename... Args>
    void warn(spdlog::format_string_t<Args...> fmt, Args&&... args)
    {
        logger_->warn(fmt, std::forward<Args>(args)...);
    }

    void log_model_information(const BayesModel& model, MCMCParams params);

    void log_result(const MCMCResult& result, const BayesModel& model);

   private:
    void log_iter_header(const BayesModel& model);
    std::shared_ptr<spdlog::logger> logger_;
};
}  // namespace detail
}  // namespace gelex
