#include "gelex/model/bayes/prior_manager.h"

#include <cmath>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include "Eigen/Core"

#include "../src/model/bayes/bayes_effects.h"
#include "gelex/error.h"
#include "gelex/logger.h"
#include "gelex/model/effects.h"

namespace gelex
{

PriorManager::PriorManager(BayesAlphabet alphabet)
    : alphabet_(alphabet),
      shared_marker_variance_(is_shared_marker_variance(alphabet))
{
}

auto PriorManager::default_prior(BayesModel& model)
    -> std::expected<void, Error>
{
    const double y_var = model.phenotype_var();

    if (auto* effect = model.additive(); effect != nullptr)
    {
        auto default_pi = default_mixture_prop(alphabet_);
        effect->pi = default_pi;
        auto result = set_variance(*effect, 0.5 * y_var);
        if (!result)
        {
            return std::unexpected(result.error());
        }
    }

    if (auto& effects = model.random(); !effects.empty())
    {
        for (auto& effect : effects)
        {
            auto result = set_variance(effect, 0.05 * y_var);
            if (!result)
            {
                return std::unexpected(result.error());
            }
        }
    }

    auto result = set_variance(model.residual(), 0.5 * y_var);
    if (!result)
    {
        return std::unexpected(result.error());
    }

    return {};
}

double calculate_lambda(double nu)
{
    const double pi = M_PI;
    const double term1 = 2.0 * std::sqrt((nu - 2.0) / pi);
    const double gamma_ratio
        = std::tgamma((nu + 1.0) / 2.0) / std::tgamma(nu / 2.0);
    const double term2 = gamma_ratio / (nu - 1.0);
    return term1 * term2;
}

auto PriorManager::set_dominant_ratio_prior(
    BayesModel& model,
    double mu,
    double variance) -> std::expected<void, Error>
{
    auto* dom_eff = model.dominant();

    const auto& freq_p_2 = bayes::get_means(model.additive()->design_matrix);

    dom_eff->wj = (1 - freq_p_2.array()).matrix();
    dom_eff->ratio_mean = mu;
    dom_eff->ratio_variance = variance;
    return {};
}

auto PriorManager::set_dominant_ratio_prior(
    BayesModel& model,
    double nu,
    double h2,
    double d_by_a,
    double i_prop) -> std::expected<void, Error>
{
    auto var_prior_result = init_variance_prop(model, d_by_a, h2, i_prop);
    if (!var_prior_result)
    {
        return std::unexpected(var_prior_result.error());
    }
    const auto [Va, Vd, I] = var_prior_result.value();

    const auto& add_mat
        = bayes::get_matrix_ref(model.additive()->design_matrix);
    const auto& dom_mat
        = bayes::get_matrix_ref(model.dominant()->design_matrix);
    const auto M = static_cast<double>(add_mat.cols());
    const Eigen::VectorXd p_freq
        = bayes::get_means(model.additive()->design_matrix).array() / 2.0;
    const Eigen::VectorXd q_freq = 1.0 - p_freq.array();
    const auto H_obs = bayes::get_means(model.dominant()->design_matrix);

    const double H_bar = H_obs.mean();
    const double H2_bar = H_obs.array().square().mean();

    const double gamma_M_prime
        = (H_obs.array() * (q_freq - p_freq).array().square()).mean();

    auto lambda = calculate_lambda(nu);

    if (std::abs(I) < 1e-9 || std::abs(H2_bar) < 1e-9)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::InvalidArgument,
                .message = "Inbreeding depression or H2_bar is zero, division "
                           "by zero."});
    }
    const double ratio
        = ((Vd * M * H_bar * H_bar * lambda * lambda) / (I * I * H2_bar)) - 1.0;

    if (ratio < 0)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::InvalidArgument,
                .message
                = "Calculated ratio (sigma_delta^2/mu_delta^2) is negative. "
                  "Check input Va, Vd, I."});
    }

    if (std::abs(H2_bar) < 1e-9 || std::abs(lambda) < 1e-9 || M == 0
        || std::abs(H_bar) < 1e-9)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::InvalidArgument,
                .message
                = "Division by zero when calculating quadratic coefficients."});
    }
    const double a_coeff = Va - ((gamma_M_prime * Vd) / H2_bar);

    // const double b_coeff = 0.0;
    const double c_coeff = -(I * I) / (lambda * lambda * M * H_bar);

    const double discriminant = -4.0 * a_coeff * c_coeff;  // b_coeff^2 项消失
    if (discriminant < 0)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::InvalidArgument,
                .message
                = "Discriminant is negative. No real solution for mu_delta."});
    }

    if (std::abs(a_coeff) < 1e-9)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::InvalidArgument,
                .message = "Quadratic coefficient 'a' is zero, cannot solve."});
    }

    const double mu_delta = std::sqrt(discriminant) / (2.0 * a_coeff);
    if (mu_delta < 0)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::InvalidArgument,
                .message
                = "Calculated mu_delta is negative, review inputs or model "
                  "assumptions."});
    }

    const double sigma_delta_sq = mu_delta * mu_delta * ratio;

    const double denominator_s2
        = (sigma_delta_sq + mu_delta * mu_delta) * M * H2_bar;
    if (std::abs(denominator_s2) < 1e-9)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::InvalidArgument,
                .message = "Division by zero when calculating s^2."});
    }
    const double s2 = (Vd / denominator_s2) * (nu - 2.0) / nu;

    auto* dom_effect = model.dominant();
    auto* add_effect = model.additive();

    add_effect->prior = {nu, s2};
    add_effect->init_marker_variance = s2;
    dom_effect->ratio_mean = mu_delta;
    dom_effect->ratio_variance = sigma_delta_sq;
    auto logger = logging::get();
    logger->info(
        "Set dominant effect prior: nu = {}, s^2 = {:.6f}, mu_delta = "
        "{:.6f}, sigma_delta^2 = {:.6f}",
        nu,
        s2,
        mu_delta,
        sigma_delta_sq);
    return {};
}

auto PriorManager::init_variance_prop(
    const BayesModel& model,
    double h2,
    double d_by_a,
    double i_prop) -> std::expected<std::tuple<double, double, double>, Error>
{
    double additive = 0.0;
    double dominance = 0.0;
    double inbreeding = 0.0;

    if (h2 < 0.0 || h2 > 1.0)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidArgument,
                fmt::format(
                    "Heritability must be between 0 and 1, got {}", h2)});
    }
    if (i_prop < 0.0 || i_prop > 1.0)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidArgument,
                fmt::format(
                    "Interaction proportion must be between 0 and 1, got {}",
                    i_prop)});
    }

    const double total_var = model.phenotype_var();
    additive = h2 * total_var;
    dominance = additive * d_by_a;
    inbreeding = sqrt(total_var) * i_prop;

    return std::make_tuple(additive, dominance, inbreeding);
}

auto PriorManager::set_mixture_prop(
    BayesModel& model,
    std::span<double> mixture_prop) -> std::expected<void, Error>
{
    Eigen::VectorXd prop = Eigen::Map<const Eigen::VectorXd>(
        mixture_prop.data(), mixture_prop.size());

    if (std::abs(prop.sum() - 1.0) > 1e-8)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidArgument,
                "Mixture proportions must sum to 1"});
    }

    if (!is_mixture_model(alphabet_))
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidArgument,
                "non-mixture model does not support mixture proportions"});
    }

    int expected_components = get_mixture_components(alphabet_);
    if (prop.size() != expected_components)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidArgument,
                fmt::format(
                    "Mixture proportion size must be {} for {}, got {}",
                    expected_components,
                    alphabet_,
                    prop.size())});
    }

    model.additive()->pi = prop;
    return {};
}

Eigen::VectorXd PriorManager::default_mixture_prop(BayesAlphabet type)
{
    switch (type)
    {
        case BayesAlphabet::A:
        case BayesAlphabet::RR:
            return Eigen::VectorXd{{0.0, 1.0}};

        case BayesAlphabet::B:
        case BayesAlphabet::C:
        case BayesAlphabet::Bpi:
        case BayesAlphabet::Cpi:
            return Eigen::VectorXd{{0.95, 0.05}};
        case BayesAlphabet::R:
            return Eigen::VectorXd{{0.95, 0.05, 0.02, 0.02, 0.01}};
    }
}

bool PriorManager::is_shared_marker_variance(BayesAlphabet type)
{
    switch (type)
    {
        case BayesAlphabet::RR:
        case BayesAlphabet::C:
        case BayesAlphabet::Cpi:
        case BayesAlphabet::R:
            return true;

        case BayesAlphabet::A:
        case BayesAlphabet::B:
        case BayesAlphabet::Bpi:
            return false;
    }
}
}  // namespace gelex
