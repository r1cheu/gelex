#pragma once

#include <random>

#include "../src/model/bayes/samplers/additive.h"
#include "../src/model/bayes/samplers/common.h"
#include "gelex/model/bayes/model.h"
#include "model/bayes/samplers/dominant.h"
#include "model/bayes/samplers/pi.h"

namespace gelex
{

template <typename T>
concept Sampler = requires(
    const T& op,
    const BayesModel& model,
    BayesState& status,
    std::mt19937_64& rng) {
    { op(model, status, rng) } -> std::same_as<void>;
};

template <Sampler... Samplers>
class TraitModel
{
   public:
    auto operator()(
        const BayesModel& model,
        BayesState& status,
        std::mt19937_64& rng) const -> void
    {
        std::apply(
            [&](const auto&... sampler) { (sampler(model, status, rng), ...); },
            samplers_);
    }

   private:
    std::tuple<Samplers...> samplers_{};
};

template <Sampler... Samplers>
using TraitBasicDefault = TraitModel<
    detail::CommonSampler::Fixed,
    detail::CommonSampler::Random,
    Samplers...,
    detail::CommonSampler::Residual>;

using BayesRR = TraitBasicDefault<detail::AdditiveSampler::RR>;
using BayesA = TraitBasicDefault<detail::AdditiveSampler::A>;
using BayesB = TraitBasicDefault<detail::AdditiveSampler::B>;
using BayesC = TraitBasicDefault<detail::AdditiveSampler::C>;
using BayesBpi = TraitBasicDefault<detail::AdditiveSampler::B, detail::Pi::Pi>;
using BayesCpi = TraitBasicDefault<detail::AdditiveSampler::C, detail::Pi::Pi>;

using BayesRRD = TraitBasicDefault<
    detail::AdditiveSampler::RRD,
    detail::DominantSampler::Coeff>;

}  // namespace gelex
