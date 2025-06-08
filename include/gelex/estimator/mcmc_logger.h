#pragma once
#include <memory>

#include <barkeep.h>
#include <spdlog/logger.h>

#include "gelex/estimator/mcmc_params.h"
#include "gelex/estimator/mcmc_result.h"
#include "gelex/model/bayes.h"

namespace gelex
{

class MCMCLogger
{
   public:
    MCMCLogger();
    void set_verbose(bool verbose);

    void log_model_information(const BayesModel& model, MCMCParams params);
    void log_burnin_finished();

    void log_result(const MCMCResult& result);
    static std::shared_ptr<barkeep::CompositeDisplay> progress_bar(
        std::vector<size_t>& idxs,
        size_t total);

   private:
    void log_iter_header(const BayesModel& model);
    std::shared_ptr<spdlog::logger> logger_;
};

}  // namespace gelex
