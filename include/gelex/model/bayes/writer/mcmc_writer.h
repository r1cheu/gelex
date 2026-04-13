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

#ifndef GELEX_MODEL_BAYES_WRITER_MCMC_WRITER_H_
#define GELEX_MODEL_BAYES_WRITER_MCMC_WRITER_H_

#include <optional>
#include <string_view>
#include <vector>

#include <Eigen/Core>

#include "gelex/io/binary_writer.h"
#include "gelex/model/bayes/model.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class BayesModel;

namespace bayes
{
class Priors;
}

namespace mcmc
{

class Writer
{
   public:
    Writer(
        const BayesModel& model,
        const bayes::Priors& priors,
        std::string_view prefix,
        Eigen::Index n_records);
    Writer(const Writer&) = delete;
    Writer(Writer&&) noexcept = default;
    auto operator=(const Writer&) -> Writer& = delete;
    auto operator=(Writer&&) noexcept -> Writer& = default;
    ~Writer() = default;

    void write(const mcmc::State& state);
    void finalize();

   private:
    struct RandomHandles
    {
        detail::SectionHandle<double> coeffs;
        detail::SectionHandle<double> variance;
    };

    struct GeneticHandles
    {
        EffectType section_effect{};
        detail::SectionHandle<double> coeffs{};
        detail::SectionHandle<double> variance{};
        std::optional<detail::SectionHandle<int8_t>> mixture_tracker;
        std::optional<detail::SectionHandle<double>> pi;
        std::optional<detail::SectionHandle<int8_t>> sign_tracker;
        std::optional<detail::SectionHandle<double>> positive_prob;
    };

    detail::BinaryWriter writer_;
    detail::SectionHandle<double> fixed_coeffs_{};
    std::vector<RandomHandles> random_;
    std::vector<GeneticHandles> genetic_;
    detail::SectionHandle<double> residual_variance_{};
};

}  // namespace mcmc

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_WRITER_MCMC_WRITER_H_
