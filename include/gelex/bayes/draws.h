// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_DRAWS_H_
#define GELEX_BAYES_DRAWS_H_

#include <cassert>
#include <cstdint>
#include <ranges>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/bayes/genetic/factory.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/serialization_ids.h"
#include "gelex/bayes/state.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/dense_writer.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

class RandomEffectDraws
{
   public:
    RandomEffectDraws(
        DenseStream<float> coefficients,
        DenseStream<double> variance)
        : coefficients_{std::move(coefficients)}, variance_{std::move(variance)}
    {
    }

    auto append(const RandomEffectState& state) -> void
    {
        coefficients_ << state.coefficients.cast<float>().eval();
        variance_ << state.variance;
    }

   private:
    DenseStream<float> coefficients_;
    DenseStream<double> variance_;
};

GELEX_NAMESPACE_BEGIN(detail)
auto make_random_draws(
    const BayesModel& model,
    DenseWriter& writer,
    std::uint64_t draw_count) -> std::vector<RandomEffectDraws>;
GELEX_NAMESPACE_END(detail)

template <typename GeneticPrior>
class BayesDraws
{
   public:
    using genetic_draws_type = genetic_draws_t<genetic_state_t<GeneticPrior>>;

    BayesDraws(
        const BayesState<GeneticPrior>& state,
        const BayesModel& model,
        std::string_view output_path,
        std::uint64_t draw_count)
        : writer_{output_path},
          fixed_{writer_.reserve<double>(
              fixed_coefficients_id,
              BinaryShape{
                  static_cast<std::uint64_t>(state.fixed().coefficients.size()),
                  draw_count})},
          random_{detail::make_random_draws(model, writer_, draw_count)},
          genetic_{make_draws(state.genetic(), writer_, draw_count)},
          residual_{writer_.reserve<double>(
              residual_variance_id,
              BinaryShape{1, draw_count})}
    {
    }
    BayesDraws(const BayesDraws&) = delete;
    BayesDraws(BayesDraws&&) = delete;
    auto operator=(const BayesDraws&) -> BayesDraws& = delete;
    auto operator=(BayesDraws&&) -> BayesDraws& = delete;
    ~BayesDraws() noexcept = default;

    auto append(const BayesState<GeneticPrior>& state) -> void
    {
        assert(random_.size() == state.random().size());
        fixed_ << state.fixed().coefficients;
        for (auto&& [draws, random_state] :
             std::views::zip(random_, state.random()))
        {
            draws.append(random_state);
        }
        genetic_.for_each(
            [&]<GeneticMode Mode>(auto& draws)
            { draws.append(state.genetic().template get<Mode>()); });
        if constexpr (requires { genetic_.joint(); })
        {
            genetic_.joint().append(state.genetic().joint());
        }
        residual_ << state.residual().variance;
    }

    // Publishes the file; requires every reserved draw to be appended.
    auto close() -> void { writer_.close(); }

   private:
    // Every leaf DenseStream borrows this address.
    DenseWriter writer_;
    DenseStream<double> fixed_;
    std::vector<RandomEffectDraws> random_;
    genetic_draws_type genetic_;
    DenseStream<double> residual_;
};

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_DRAWS_H_
