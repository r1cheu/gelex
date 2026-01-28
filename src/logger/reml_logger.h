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
