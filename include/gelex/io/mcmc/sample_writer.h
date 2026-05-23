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

#ifndef GELEX_IO_MCMC_SAMPLE_WRITER_H_
#define GELEX_IO_MCMC_SAMPLE_WRITER_H_

#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/io/detail/binary_writer.h"
#include "gelex/model/bayes/legacy_method.h"
#include "gelex/model/bayes/model.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class BayesModel;
class BayesState;

namespace mcmc
{

class Writer
{
   public:
    using RecordHandle = std::variant<
        io::detail::SectionHandle<float>,
        io::detail::SectionHandle<double>,
        io::detail::SectionHandle<int>>;

    Writer(
        const BayesModel& model,
        const bayes::LegacyBayesMethod& method,
        std::string_view prefix,
        Eigen::Index n_records);
    Writer(
        const BayesState& state,
        std::string_view prefix,
        Eigen::Index n_records);
    Writer(const Writer&) = delete;
    Writer(Writer&&) = delete;
    auto operator=(const Writer&) -> Writer& = delete;
    auto operator=(Writer&&) -> Writer& = delete;
    ~Writer() = default;

    void write(const mcmc::State& state);
    void write(const BayesState& state);

   private:
    struct RandomHandles
    {
        io::detail::SectionHandle<double> coeffs;
        io::detail::SectionHandle<double> variance;
    };

    struct GeneticHandles
    {
        EffectType section_effect{};
        io::detail::SectionHandle<double> coeffs{};
        io::detail::SectionHandle<double> variance{};
        std::optional<io::detail::SectionHandle<int8_t>> mixture_tracker;
        std::optional<io::detail::SectionHandle<double>> pi;
        std::optional<io::detail::SectionHandle<int8_t>> sign_tracker;
        std::optional<io::detail::SectionHandle<double>> positive_prob;
    };

    io::detail::BinaryWriter writer_;
    std::unordered_map<std::string, RecordHandle> record_handles_;
    io::detail::SectionHandle<double> fixed_coeffs_{};
    std::vector<RandomHandles> random_;
    std::vector<GeneticHandles> genetic_;
    io::detail::SectionHandle<double> residual_variance_{};
};

}  // namespace mcmc

}  // namespace gelex

#endif  // GELEX_IO_MCMC_SAMPLE_WRITER_H_
