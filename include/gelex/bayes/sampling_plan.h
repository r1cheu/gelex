// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_SAMPLING_PLAN_H_
#define GELEX_BAYES_SAMPLING_PLAN_H_

#include <cstdint>

namespace gelex
{

// Iteration schedule of one MCMC chain: every draw after burn-in at the
// thinning interval is retained, so thin must divide iterations - burn_in.
class SamplingPlan
{
   public:
    SamplingPlan(int iterations, int burn_in, int thin, int seed);

    [[nodiscard]] auto iterations() const noexcept -> int
    {
        return iterations_;
    }
    [[nodiscard]] auto burn_in() const noexcept -> int { return burn_in_; }
    [[nodiscard]] auto thin() const noexcept -> int { return thin_; }
    [[nodiscard]] auto seed() const noexcept -> int { return seed_; }

    [[nodiscard]] auto retains(int iteration) const noexcept -> bool
    {
        return iteration >= burn_in_ && (iteration + 1 - burn_in_) % thin_ == 0;
    }

    [[nodiscard]] auto draw_count() const noexcept -> std::uint64_t
    {
        return static_cast<std::uint64_t>((iterations_ - burn_in_) / thin_);
    }

    auto operator==(const SamplingPlan&) const noexcept -> bool = default;

   private:
    int iterations_;
    int burn_in_;
    int thin_;
    int seed_;
};

}  // namespace gelex

#endif  // GELEX_BAYES_SAMPLING_PLAN_H_
