#pragma once

#include <memory>
#include "gelex/model/bayes/prior_strategy.h"
#include "gelex/model/effects.h"
#include "types/bayes_effects.h"

namespace gelex
{
/**
 * BayesA prior strategy
 *
 * All markers have non-zero effects with common variance
 */
class BayesAPrior : public PriorSetter
{
   public:
    ~BayesAPrior() override = default;

   private:
    auto set_additive_effect_prior(
        bayes::AdditiveEffect& effect,
        const PriorConfig& prior) -> std::expected<void, Error> override;
};

/**
 * BayesA with Dominance prior strategy
 *
 * Extends BayesA with dominant effect priors
 */
class BayesADPrior : public BayesAPrior
{
   public:
    ~BayesADPrior() override = default;

   private:
    auto set_dominant_effect_prior(
        bayes::DominantEffect& effect,
        const PriorConfig& prior) -> std::expected<void, Error> override;
};

/**
 * BayesB prior strategy
 *
 * Mixture prior with some markers having zero effects
 */
class BayesBPrior : public PriorSetter
{
   public:
    ~BayesBPrior() override = default;

   private:
    auto set_additive_effect_prior(
        bayes::AdditiveEffect& effect,
        const PriorConfig& prior) -> std::expected<void, Error> override;
};

/**
 * BayesB with Dominance prior strategy
 *
 * Extends BayesB with dominant effect priors
 */
class BayesBDPrior : public BayesBPrior
{
   public:
    ~BayesBDPrior() override = default;

   private:
    auto set_dominant_effect_prior(
        bayes::DominantEffect& effect,
        const PriorConfig& prior) -> std::expected<void, Error> override;
};

/**
 * BayesC prior strategy
 *
 * Common variance for all non-zero markers
 */
class BayesCPrior : public PriorSetter
{
   public:
    ~BayesCPrior() override = default;

   private:
    auto set_additive_effect_prior(
        bayes::AdditiveEffect& effect,
        const PriorConfig& prior) -> std::expected<void, Error> override;
};

/**
 * BayesC with Dominance prior strategy
 *
 * Extends BayesC with dominant effect priors
 */
class BayesCDPrior : public BayesCPrior
{
   public:
    ~BayesCDPrior() override = default;

   private:
    auto set_dominant_effect_prior(
        bayes::DominantEffect& effect,
        const PriorConfig& prior) -> std::expected<void, Error> override;
};

/**
 * Bayesian Ridge Regression prior strategy
 *
 * All markers have non-zero effects with common variance (equivalent to
 * RR-BLUP)
 */
class BayesRRPrior : public PriorSetter
{
   public:
    ~BayesRRPrior() override = default;

   private:
    auto set_additive_effect_prior(
        bayes::AdditiveEffect& effect,
        const PriorConfig& prior) -> std::expected<void, Error> override;
};

/**
 * Bayesian Ridge Regression with Dominance prior strategy
 *
 * Extends BayesRR with dominant effect priors
 */
class BayesRRDPrior : public BayesRRPrior
{
   public:
    ~BayesRRDPrior() override = default;

   private:
    auto set_dominant_effect_prior(
        bayes::DominantEffect& effect,
        const PriorConfig& prior) -> std::expected<void, Error> override;
};

/**
 * BayesR prior strategy
 *
 * Mixture prior with scaled variances for different effect sizes
 */
class BayesRPrior : public PriorSetter
{
   public:
    ~BayesRPrior() override = default;

   private:
    auto set_additive_effect_prior(
        bayes::AdditiveEffect& effect,
        const PriorConfig& prior) -> std::expected<void, Error> override;
};

/**
 * BayesR with Dominance prior strategy
 *
 * Extends BayesR with dominant effect priors
 */
class BayesRDPrior : public BayesRPrior
{
   public:
    ~BayesRDPrior() override = default;

   private:
    auto set_dominant_effect_prior(
        bayes::DominantEffect& effect,
        const PriorConfig& prior) -> std::expected<void, Error> override;
};

/**
 * Factory function to create appropriate prior strategy based on model type
 *
 * @param type The BayesAlphabet model type
 * @return Unique pointer to the appropriate PriorSetter instance
 */
inline auto create_prior_strategy(BayesAlphabet type)
    -> std::unique_ptr<PriorSetter>
{
    using bt = BayesAlphabet;

    switch (type)
    {
        case bt::A:
            return std::make_unique<BayesAPrior>();
        case bt::Ad:
            return std::make_unique<BayesADPrior>();
        case bt::B:
        case bt::Bpi:
            return std::make_unique<BayesBPrior>();
        case bt::Bd:
        case bt::Bdpi:
            return std::make_unique<BayesBDPrior>();
        case bt::C:
        case bt::Cpi:
            return std::make_unique<BayesCPrior>();
        case bt::Cd:
        case bt::Cdpi:
            return std::make_unique<BayesCDPrior>();
        case bt::R:
            return std::make_unique<BayesRPrior>();
        case bt::Rd:
            return std::make_unique<BayesRDPrior>();
        case bt::RR:
            return std::make_unique<BayesRRPrior>();
        case bt::RRd:
            return std::make_unique<BayesRRDPrior>();
        default:
            return nullptr;
    }
}

}  // namespace gelex
