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

#ifndef GELEX_LOGGER_REML_LOGGER_BASE_H_
#define GELEX_LOGGER_REML_LOGGER_BASE_H_

#include <cstddef>
#include <memory>

#include <spdlog/logger.h>
#include "gelex/logger.h"

namespace gelex
{
class FreqModel;
class FreqState;

namespace detail
{

class RemlLoggerBase
{
   public:
    RemlLoggerBase() : logger_{gelex::logging::get()} {}
    RemlLoggerBase(const RemlLoggerBase&) = default;
    RemlLoggerBase(RemlLoggerBase&&) = delete;
    RemlLoggerBase& operator=(const RemlLoggerBase&) = default;
    RemlLoggerBase& operator=(RemlLoggerBase&&) = delete;
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
