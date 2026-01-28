/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef GELEX_LOGGER_BAYES_LOGGER_H_
#define GELEX_LOGGER_BAYES_LOGGER_H_
#include <memory>

#include <barkeep.h>
#include <spdlog/logger.h>
#include <Eigen/Core>

#include "gelex/model/bayes/model.h"
#include "gelex/types/mcmc_results.h"

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

    void log_model_information(const BayesModel& model);

    void log_result(
        const MCMCResult& result,
        const BayesModel& model,
        double elapsed_time,
        Eigen::Index samples_collected);

   private:
    void log_iter_header(const BayesModel& model);
    std::shared_ptr<spdlog::logger> logger_;
};
}  // namespace detail
}  // namespace gelex

#endif  // GELEX_LOGGER_BAYES_LOGGER_H_
