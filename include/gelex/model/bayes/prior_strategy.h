#pragma once

#include <Eigen/Core>
#include <vector>

#include "gelex/model/bayes/model.h"
#include "types/bayes_effects.h"

namespace gelex
{

/**
 * Configuration for a single effect type prior
 */
struct Prior
{
    Eigen::VectorXd
        mixture_proportions;  ///< Mixture proportions for effect sizes
    Eigen::VectorXd
        mixture_scales;   ///< Scale parameters for mixture components
    double heritability;  ///< Heritability proportion for this effect type
};

/**
 * Complete prior configuration for Bayesian models
 *
 * Contains prior settings for all effect types and variance components
 */
struct PriorConfig
{
    double phenotype_variance{0.0};  ///< Total phenotypic variance
    Prior additive{
        Eigen::VectorXd::Zero(2),
        Eigen::VectorXd(5),
        0.5};  ///< Additive effect prior configuration
    Prior dominant{
        Eigen::VectorXd::Zero(2),
        Eigen::VectorXd(5),
        0.2};  ///< Dominant effect prior configuration
    double random_variance_proportion{
        0.1};  ///< Proportion of variance for random effects
    double residual_variance_proportion{
        0.3};  ///< Proportion of variance for residuals
};

/**
 * Abstract base class for setting priors in Bayesian models
 *
 * Implements the Strategy pattern for different prior configurations
 */
class PriorSetter
{
   public:
    PriorSetter() = default;
    PriorSetter(const PriorSetter&) = default;
    PriorSetter(PriorSetter&&) = default;
    PriorSetter& operator=(const PriorSetter&) = default;
    PriorSetter& operator=(PriorSetter&&) = default;
    virtual ~PriorSetter() = default;

    virtual auto operator()(BayesModel& model, const PriorConfig& config)
        -> std::expected<void, Error>;

   protected:
    static auto set_random_effect_prior(
        double variance,
        bayes::RandomEffect& effect) -> void;

    static auto compute_init_marker_variance(
        double target_variance,
        const Eigen::Ref<const Eigen::MatrixXd>& design_matrix,
        double non_zero_marker_proption) -> std::expected<double, Error>;
    // use for unstandardized genotype matrix
    static auto compute_init_marker_variance(
        double target_variance,
        const Eigen::Ref<const Eigen::VectorXd>& genetic_variance,
        double non_zero_marker_proption) -> std::expected<double, Error>;

   private:
    virtual auto set_additive_effect_prior(
        bayes::AdditiveEffect& effect,
        const PriorConfig& prior) -> std::expected<void, Error>
        = 0;
    virtual auto set_dominant_effect_prior(
        bayes::DominantEffect& effect,
        const PriorConfig& prior) -> std::expected<void, Error>
    {
        return {};
    };
};
}  // namespace gelex
