// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_DRAWS_H_
#define GELEX_BAYES_DRAWS_H_

#include <Eigen/Core>
#include <cassert>
#include <cstdint>
#include <fmt/format.h>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/bayes/genetic/draw_traits.h"
#include "gelex/bayes/genetic/factory.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/serialization_ids.h"
#include "gelex/bayes/state.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/csc_writer.h"
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

    auto operator<<(const RandomEffectState& state) -> RandomEffectDraws&
    {
        scratch_ = state.coefficients.cast<float>();
        coefficients_ << scratch_;
        variance_ << state.variance;
        return *this;
    }

   private:
    DenseStream<float> coefficients_;
    DenseStream<double> variance_;
    // Level-length float conversion buffer, reused across draws.
    Eigen::VectorXf scratch_;
};

// Marker-level sparse draws live next to the dense file.
inline auto sparse_draws_path(std::string_view output_path) -> std::string
{
    return fmt::format("{}.csc", output_path);
}

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

    // Dense payloads go to `output_path`, marker-level sparse payloads to
    // sparse_draws_path(output_path).
    BayesDraws(
        const BayesState<GeneticPrior>& state,
        const BayesModel& model,
        std::string_view output_path,
        std::uint64_t draw_count)
        : dense_{output_path},
          sparse_{sparse_draws_path(output_path)},
          fixed_{dense_.reserve<double>(
              fixed_coefficients_id,
              BinaryShape{
                  static_cast<std::uint64_t>(state.fixed().coefficients.size()),
                  draw_count})},
          random_{detail::make_random_draws(model, dense_, draw_count)},
          genetic_{make_draws(
              state.genetic(),
              DrawWriters{.dense = dense_, .sparse = sparse_},
              draw_count)},
          residual_{dense_.reserve<double>(
              residual_variance_id,
              BinaryShape{1, draw_count})}
    {
    }
    BayesDraws(const BayesDraws&) = delete;
    BayesDraws(BayesDraws&&) = delete;
    auto operator=(const BayesDraws&) -> BayesDraws& = delete;
    auto operator=(BayesDraws&&) -> BayesDraws& = delete;
    ~BayesDraws() noexcept = default;

    auto operator<<(const BayesState<GeneticPrior>& state) -> BayesDraws&
    {
        assert(random_.size() == state.random().size());
        fixed_ << state.fixed().coefficients;
        for (auto&& [draws, random_state] :
             std::views::zip(random_, state.random()))
        {
            draws << random_state;
        }
        genetic_.for_each([&]<GeneticMode Mode>(auto& draws)
                          { draws << state.genetic().template get<Mode>(); });
        if constexpr (requires { genetic_.joint(); })
        {
            genetic_.joint() << state.genetic().joint();
        }
        residual_ << state.residual().variance;
        return *this;
    }

    // Publishes both files; requires every reserved draw to be appended.
    // The dense file carries every model, so a short run fails there before
    // anything is published.
    auto close() -> void
    {
        dense_.close();
        sparse_.close();
    }

   private:
    // Every leaf stream borrows one of these addresses.
    DenseWriter dense_;
    CscWriter sparse_;
    DenseStream<double> fixed_;
    std::vector<RandomEffectDraws> random_;
    genetic_draws_type genetic_;
    DenseStream<double> residual_;
};

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_DRAWS_H_
