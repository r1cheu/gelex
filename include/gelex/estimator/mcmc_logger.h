#pragma once
#include <memory>

#include <spdlog/logger.h>

#include "gelex/estimator/mcmc_params.h"
#include "gelex/model/bayes.h"

namespace gelex
{

class MCMCLogger
{
   public:
    MCMCLogger();
    void set_verbose(bool verbose);

    void log_model_information(const Bayes& model, MCMCParams params);
    void
    log_iteration(size_t iter, const Bayes& model, std::string_view duartion);
    void log_burnin_finished();

   private:
    std::shared_ptr<spdlog::logger> logger_;
};

}  // namespace gelex
