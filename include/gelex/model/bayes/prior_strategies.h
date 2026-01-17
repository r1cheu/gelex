#ifndef GELEX_MODEL_BAYES_PRIOR_STRATEGIES_H_
#define GELEX_MODEL_BAYES_PRIOR_STRATEGIES_H_

#include <memory>
#include "../src/types/bayes_effects.h"
#include "gelex/model/bayes/prior_strategy.h"
#include "gelex/model/effects.h"

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
        const PriorConfig& prior) -> void override;
};

/**
 * BayesA with Dominance prior strategy
 *
 * Extends BayesA with dominant effect priors
 */
class BayesAdPrior : public BayesAPrior
{
   public:
    ~BayesAdPrior() override = default;

   private:
    auto set_dominant_effect_prior(
        bayes::DominantEffect& effect,
        const PriorConfig& prior) -> void override;
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
        const PriorConfig& prior) -> void override;
};

/**
 * BayesBpi prior strategy
 *
 * Mixture prior with some markers having zero effects and estimating pi
 */
class BayesBpiPrior : public BayesBPrior
{
   public:
    ~BayesBpiPrior() override = default;

   private:
    auto set_additive_effect_prior(
        bayes::AdditiveEffect& effect,
        const PriorConfig& prior) -> void override;
};

/**
 * BayesB with Dominance prior strategy
 *
 * Extends BayesB with dominant effect priors
 */
class BayesBdPrior : public BayesBPrior
{
   public:
    ~BayesBdPrior() override = default;

   private:
    auto set_dominant_effect_prior(
        bayes::DominantEffect& effect,
        const PriorConfig& prior) -> void override;
};

/**
 * BayesBdpi prior strategy
 *
 * Extends BayesB with dominant effect priors and estimating pi for both
 * additive and dominant effects
 */
class BayesBdpiPrior : public BayesBpiPrior
{
   public:
    ~BayesBdpiPrior() override = default;

   private:
    auto set_dominant_effect_prior(
        bayes::DominantEffect& effect,
        const PriorConfig& prior) -> void override;
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
        const PriorConfig& prior) -> void override;
};

/**
 * BayesCpi prior strategy
 *
 * Common variance for all non-zero markers and estimating pi
 */
class BayesCpiPrior : public BayesCPrior
{
   public:
    ~BayesCpiPrior() override = default;

   private:
    auto set_additive_effect_prior(
        bayes::AdditiveEffect& effect,
        const PriorConfig& prior) -> void override;
};

/**
 * BayesC with Dominance prior strategy
 *
 * Extends BayesC with dominant effect priors
 */
class BayesCdPrior : public BayesCPrior
{
   public:
    ~BayesCdPrior() override = default;

   private:
    auto set_dominant_effect_prior(
        bayes::DominantEffect& effect,
        const PriorConfig& prior) -> void override;
};

/**
 * BayesCdpi prior strategy
 *
 * Extends BayesC with dominant effect priors and estimating pi for both
 * additive and dominant effects
 */
class BayesCdpiPrior : public BayesCpiPrior
{
   public:
    ~BayesCdpiPrior() override = default;

   private:
    auto set_dominant_effect_prior(
        bayes::DominantEffect& effect,
        const PriorConfig& prior) -> void override;
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
        const PriorConfig& prior) -> void override;
};

/**
 * Bayesian Ridge Regression with Dominance prior strategy
 *
 * Extends BayesRR with dominant effect priors
 */
class BayesRRdPrior : public BayesRRPrior
{
   public:
    ~BayesRRdPrior() override = default;

   private:
    auto set_dominant_effect_prior(
        bayes::DominantEffect& effect,
        const PriorConfig& prior) -> void override;
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
        const PriorConfig& prior) -> void override;
};

/**
 * BayesR with Dominance prior strategy
 *
 * Extends BayesR with dominant effect priors
 */
class BayesRdPrior : public BayesRPrior
{
   public:
    ~BayesRdPrior() override = default;

   private:
    auto set_dominant_effect_prior(
        bayes::DominantEffect& effect,
        const PriorConfig& prior) -> void override;
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
            return std::make_unique<BayesAdPrior>();
        case bt::B:
            return std::make_unique<BayesBPrior>();
        case bt::Bpi:
            return std::make_unique<BayesBpiPrior>();
        case bt::Bd:
            return std::make_unique<BayesBdPrior>();
        case bt::Bdpi:
            return std::make_unique<BayesBdpiPrior>();
        case bt::C:
            return std::make_unique<BayesCPrior>();
        case bt::Cpi:
            return std::make_unique<BayesCpiPrior>();
        case bt::Cd:
            return std::make_unique<BayesCdPrior>();
        case bt::Cdpi:
            return std::make_unique<BayesCdpiPrior>();
        case bt::R:
            return std::make_unique<BayesRPrior>();
        case bt::Rd:
            return std::make_unique<BayesRdPrior>();
        case bt::RR:
            return std::make_unique<BayesRRPrior>();
        case bt::RRd:
            return std::make_unique<BayesRRdPrior>();
        default:
            return nullptr;
    }
}

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_PRIOR_STRATEGIES_H_
