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

#ifndef GELEX_CLI_FIT_REPORTER_H_
#define GELEX_CLI_FIT_REPORTER_H_

#include <string_view>

#include <Eigen/Core>

#include "gelex/infra/logging/fit_event.h"

namespace gelex
{
struct PosteriorSummary;

namespace stats::detail
{
struct ScaledInvChiSqParams;
}  // namespace stats::detail

namespace bayes
{
class GeneticPrior;
class BayesPrior;
class VarianceParameter;
}  // namespace bayes

enum class GeneticMode : uint8_t;
}  // namespace gelex

namespace gelex::cli
{

class FitReporter
{
   public:
    auto on_event(const FitPriorSetEvent& event) const -> void;
    auto on_event(const FitResultsSavedEvent& event) const -> void;

   protected:
    FitReporter() = default;

    static auto print_random_prior(const bayes::VarianceParameter& parameter)
        -> void;
    static auto print_genetic_prior(const bayes::GeneticPrior& prior) -> void;
    static auto print_residual_prior(const bayes::VarianceParameter& parameter)
        -> void;
    static auto print_variance_prior(
        const stats::detail::ScaledInvChiSqParams& prior,
        double init_variance) -> void;
    static auto print_summary_row(
        std::string_view name,
        const PosteriorSummary& summary,
        Eigen::Index index = 0) -> void;
};

}  // namespace gelex::cli

#endif  // GELEX_CLI_FIT_REPORTER_H_
