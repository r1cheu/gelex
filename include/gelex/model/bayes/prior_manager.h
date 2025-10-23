#pragma once

#include <Eigen/Core>
#include <tuple>

#include "gelex/model/bayes/model.h"
#include "gelex/model/effects.h"

namespace gelex
{

class PriorManager
{
   public:
    explicit PriorManager(BayesAlphabet alphabet);

    auto default_prior(BayesModel& model) -> std::expected<void, Error>;

    template <typename Effect>
    auto set_variance(Effect& effect, double target_variance)
        -> std::expected<void, Error>
    {
        if constexpr (std::is_same_v<Effect, bayes::AdditiveEffect>)
        {
            const double genetic_var = std::visit(
                [](const auto& s) { return s.variance().sum(); },
                effect.design_matrix);

            const double init_marker_variance
                = target_variance / genetic_var / (1 - effect.pi(0));

            effect.init_marker_variance = init_marker_variance;
            effect.prior = {4.0, 0.5 * init_marker_variance};

            effect.marker_variance_size
                = shared_marker_variance_
                      ? 1
                      : bayes::get_cols(effect.design_matrix);
        }
        else if constexpr (std::is_same_v<
                               Effect,
                               std::vector<bayes::RandomEffect>>)
        {
            const size_t num_effects = effect.size();
            const double init_effect_variance
                = target_variance / static_cast<double>(num_effects);
            for (auto& eff : effect)
            {
                eff.prior = {4.0, 0.5 * init_effect_variance};
                eff.init_variance = init_effect_variance;
            }
        }
        else if constexpr (std::is_same_v<Effect, bayes::Residual>)
        {
            effect.prior = {4, 0};
            effect.init_variance = target_variance;
        }
        else
        {
            return std::unexpected(
                Error{
                    ErrorCode::InvalidArgument,
                    "Unsupported effect type for variance prior"});
        }
        return {};
    }
    auto set_dominant_ratio_prior(BayesModel& model, double mu, double variance)
        -> std::expected<void, Error>;
    auto set_dominant_ratio_prior(
        BayesModel& model,
        double nu,
        double h2,
        double d_by_a,
        double i_prop) -> std::expected<void, Error>;
    auto set_mixture_prop(BayesModel& model, std::span<double> mixture_prop)
        -> std::expected<void, Error>;

   private:
    BayesAlphabet alphabet_;
    bool shared_marker_variance_{false};

    auto init_variance_prop(
        const BayesModel& model,
        double h2,
        double d_by_a,
        double i_prop)
        -> std::expected<std::tuple<double, double, double>, Error>;

    static Eigen::VectorXd default_mixture_prop(BayesAlphabet type);
    static bool is_shared_marker_variance(BayesAlphabet type);
};

}  // namespace gelex
