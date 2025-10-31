#include "gelex/model/bayes/prior_strategy.h"

#include <expected>
#include "gelex/error.h"
#include "gelex/model/bayes/prior_constants.h"

namespace gelex
{

auto PriorSetter::operator()(BayesModel& model, const PriorConfig& config)
    -> std::expected<void, Error>
{
    // Set additive effect prior
    if (model.additive() != nullptr)
    {
        auto additive_result
            = set_additive_effect_prior(*model.additive(), config);
        if (!additive_result)
        {
            return std::unexpected(additive_result.error());
        }
    }

    // Set dominant effect prior
    if (model.dominant() != nullptr)
    {
        auto dominant_result
            = set_dominant_effect_prior(*model.dominant(), config);
        if (!dominant_result)
        {
            return std::unexpected(dominant_result.error());
        }
    }

    const auto num_random_effect = static_cast<double>(model.random().size());
    if (!model.random().empty())
    {
        const double target_variance
            = config.random_variance_proportion * config.phenotype_variance;
        const double per_effect_variance = target_variance / num_random_effect;

        for (auto& effect : model.random())
        {
            set_random_effect_prior(per_effect_variance, effect);
        }
    }

    // Set residual prior
    {
        const double target_variance
            = config.residual_variance_proportion * config.phenotype_variance;
        model.residual().prior
            = {prior_constants::RESIDUAL_VARIANCE_SHAPE,
               prior_constants::RESIDUAL_VARIANCE_SCALE};
        model.residual().init_variance = target_variance;
    }

    return {};
}

auto PriorSetter::set_random_effect_prior(
    double variance,
    bayes::RandomEffect& effect) -> void
{
    effect.prior
        = {prior_constants::RANDOM_EFFECTS_SHAPE,
           prior_constants::RANDOM_EFFECTS_SCALE};
    effect.init_variance = variance;
}

auto PriorSetter::compute_init_marker_variance(
    double target_variance,
    const Eigen::Ref<const Eigen::MatrixXd>& design_matrix,
    double non_zero_marker_proption) -> std::expected<double, Error>
{
    if (target_variance <= 0.0)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::InvalidArgument,
                .message = "Target variance must be positive"});
    }
    if (non_zero_marker_proption <= 0.0)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::InvalidArgument,
                .message = "Non-zero marker proportion must be positive"});
    }

    auto num_snps = static_cast<double>(design_matrix.cols());
    auto num_non_zero_snps = num_snps * non_zero_marker_proption;

    if (num_non_zero_snps <= 0.0)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::InvalidArgument,
                .message = "Number of non-zero SNPs must be positive"});
    }

    return target_variance / num_non_zero_snps;
}

auto PriorSetter::compute_init_marker_variance(
    double target_variance,
    const Eigen::Ref<const Eigen::VectorXd>& genetic_variance,
    double non_zero_marker_proption) -> std::expected<double, Error>
{
    if (target_variance <= 0.0)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::InvalidArgument,
                .message = "Target variance must be positive"});
    }
    if (non_zero_marker_proption <= 0.0)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::InvalidArgument,
                .message = "Non-zero marker proportion must be positive"});
    }

    auto total_genetic_variance = genetic_variance.sum();
    auto non_zero_marker_variance
        = total_genetic_variance * non_zero_marker_proption;

    if (non_zero_marker_variance <= 0.0)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::InvalidArgument,
                .message = "Non-zero marker variance must be positive"});
    }

    return target_variance / non_zero_marker_variance;
}
}  // namespace gelex
