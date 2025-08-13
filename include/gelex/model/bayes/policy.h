#pragma once
#include <fmt/format.h>
#include <armadillo>
#include <cassert>
#include <vector>

#include "gelex/dist.h"
#include "gelex/estimator/bayes/base.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/utils/formatter.h"

namespace gelex
{

template <BayesAlphabet>
struct GeneticTrait;

template <>
struct GeneticTrait<BayesAlphabet::A>
{
    static arma::dvec sigma(const arma::dmat& X)
    {
        return arma::dvec(X.n_cols, arma::fill::zeros);
    }
    static const bool estimate_pi{false};
    static arma::dvec pi() { return arma::dvec{0.0, 1.0}; }
    static std::vector<std::string>
    prior_str(double nu, double s2, const arma::dvec& pi)
    {
        std::vector<std::string> prior_strings;
        prior_strings.reserve(3);
        prior_strings.emplace_back("BayesA");
        prior_strings.emplace_back("      ├─ αᵢ ~ N(0, σ²ᵢ)");
        prior_strings.emplace_back(
            fmt::format("      └─ {}", sigma_prior("ᵢ", nu, s2)));
        return prior_strings;
    }
    static void sample(
        const bayes::GeneticEffect& design,
        bayes::GeneticEffectState& state,
        double* y_adj,
        arma::uvec& snp_tracker,
        double sigma_e,
        std::mt19937_64& rng)
    {
        arma::dvec& coeff = state.coeff;
        double* u = state.u.memptr();

        const arma::dmat& design_matrix = design.design_matrix;
        arma::dvec& sigma = state.sigma;
        ScaledInvChiSq chi_squared{design.prior};
        std::normal_distribution<double> normal{0, 1};
        const int n = static_cast<int>(design_matrix.n_rows);

        const auto norm = static_cast<double>(design_matrix.n_rows) - 1.0;

        for (size_t i = 0; i < coeff.n_elem; ++i)
        {
            const double old_i = coeff.at(i);
            chi_squared.update(old_i * old_i, 1);
            const double new_sigma = chi_squared.sample(rng);
            const double* col_i = design_matrix.colptr(i);
            const double rhs = ddot_ptr(n, col_i, y_adj) + (norm * old_i);
            const double inv_scaler = 1 / (norm + sigma_e / new_sigma);
            const double new_i = (normal(rng) * sqrt(sigma_e * inv_scaler))
                                 + (rhs * inv_scaler);

            coeff.at(i) = new_i;
            sigma.at(i) = new_sigma;

            const double diff = old_i - new_i;
            daxpy_ptr(n, diff, col_i, y_adj);
            daxpy_ptr(n, -diff, col_i, u);
        }
    }
};

template <>
struct GeneticTrait<BayesAlphabet::RR>
{
    static arma::dvec sigma(const arma::dmat& X) { return arma::dvec{0}; }

    static const bool estimate_pi{false};
    static arma::dvec pi() { return arma::dvec{0.0, 1.0}; }
    static std::vector<std::string>
    prior_str(double nu, double s2, const arma::dvec& pi)
    {
        std::vector<std::string> prior_strings;
        prior_strings.reserve(3);
        prior_strings.emplace_back("BayesRR");
        prior_strings.emplace_back("      ├─ αᵢ ~ N(0, σ²)");
        prior_strings.emplace_back(
            fmt::format("      └─ {}", sigma_prior("", nu, s2)));
        return prior_strings;
    }
    static void sample(
        const bayes::GeneticEffect& design,
        bayes::GeneticEffectState& state,
        double* y_adj,
        arma::uvec& snp_tracker,
        double sigma_e,
        std::mt19937_64& rng)
    {
        arma::dvec& coeff = state.coeff;
        double* u = state.u.memptr();
        const arma::dmat& design_matrix = design.design_matrix;
        const double sigma_g = arma::as_scalar(state.sigma);
        const double sigma_e_sqrt = sqrt(sigma_e);
        const double inv_scaler_base = sigma_e / sigma_g;
        ScaledInvChiSq chi_squared{design.prior};

        std::normal_distribution<double> normal{0, 1};
        const int n = static_cast<int>(design_matrix.n_rows);
        const auto norm = static_cast<double>(design_matrix.n_rows) - 1.0;

        for (size_t idx = 0; idx < coeff.n_elem; ++idx)
        {
            const double old_i = coeff.at(idx);
            const double inv_scaler = 1.0 / (norm + inv_scaler_base);
            const double* col_i = design_matrix.colptr(idx);
            double rhs = ddot_ptr(n, col_i, y_adj) + (norm * old_i);
            const double new_i = (normal(rng) * sigma_e_sqrt * sqrt(inv_scaler))
                                 + (rhs * inv_scaler);

            coeff.at(idx) = new_i;
            const double diff = old_i - new_i;
            daxpy_ptr(n, diff, col_i, y_adj);
            daxpy_ptr(n, -diff, col_i, u);
        }
        const double ssq = arma::dot(state.coeff, state.coeff);
        chi_squared.update(ssq, coeff.n_elem);
        state.sigma.at(0) = chi_squared.sample(rng);
    }
};

template <>
struct GeneticTrait<BayesAlphabet::B>
{
    static arma::dvec sigma(const arma::dmat& X)
    {
        return arma::zeros<arma::dvec>(X.n_cols);
    }

    static const bool estimate_pi{false};
    static arma::dvec pi() { return arma::dvec{0.95, 0.05}; }
    static std::vector<std::string>
    prior_str(double nu, double s2, const arma::dvec& pi)
    {
        std::vector<std::string> prior_strings;
        prior_strings.reserve(3);
        prior_strings.emplace_back("BayesB");
        prior_strings.emplace_back("      ├─ αᵢ ~ 0.05 N(0, σ²ᵢ) + 0.95δ₀");
        prior_strings.emplace_back(
            fmt::format("      └─ {}", sigma_prior("ᵢ", nu, s2)));
        return prior_strings;
    }
    static void sample(
        const bayes::GeneticEffect& design,
        bayes::GeneticEffectState& state,
        double* y_adj,
        arma::uvec& snp_tracker,
        double sigma_e,
        std::mt19937_64& rng)
    {
        arma::dvec logpi = arma::log(state.pi.prop);

        arma::dvec& coeff = state.coeff;
        double* u = state.u.memptr();
        const arma::dmat& design_matrix = design.design_matrix;
        arma::dvec& sigma = state.sigma;

        ScaledInvChiSq chi_squared{design.prior};
        std::normal_distribution<double> normal{0, 1};
        std::uniform_real_distribution<double> uniform{0, 1};

        const int n = static_cast<int>(design_matrix.n_rows);
        const auto norm = static_cast<double>(design_matrix.n_rows) - 1.0;

        for (size_t i = 0; i < coeff.n_elem; ++i)
        {
            const double old_i = coeff.at(i);
            chi_squared.update(old_i * old_i, 1);
            const double new_sigma = chi_squared.sample(rng);
            const double* col_i = design_matrix.colptr(i);
            const double inv_scaler = 1 / (norm + sigma_e / new_sigma);

            double rhs = ddot_ptr(n, col_i, y_adj);
            if (old_i != 0)
            {
                rhs += norm * old_i;
            }
            double logdetV = log((new_sigma * norm / sigma_e) + 1);
            double uhat = rhs * inv_scaler;

            const double l_diff = (-0.5 * (logdetV - uhat * rhs / sigma_e))
                                  + logpi.at(1) - logpi.at(0);
            const double accept_prob = 1 / (1 + exp(l_diff));
            const bool has_effect = (uniform(rng) >= accept_prob);
            snp_tracker.at(i) = has_effect ? 1 : 0;

            double new_i{};
            if (has_effect)
            {
                new_i = normal(rng) * sqrt(sigma_e * inv_scaler) + uhat;
                const double diff = old_i - new_i;
                daxpy_ptr(n, diff, col_i, y_adj);
                daxpy_ptr(n, -diff, col_i, u);
            }
            else if (old_i != 0.0)
            {
                daxpy_ptr(n, old_i, col_i, y_adj);
                daxpy_ptr(n, -old_i, col_i, u);
            }
            coeff.at(i) = new_i;
            sigma.at(i) = new_sigma;
        }
        state.pi.count.at(1) = arma::sum(snp_tracker);
        state.pi.count.at(0) = coeff.n_elem - state.pi.count.at(1);
    }
};

template <>
struct GeneticTrait<BayesAlphabet::Bpi>
{
    static arma::dvec sigma(const arma::dmat& X)
    {
        return arma::zeros<arma::dvec>(X.n_cols);
    }

    static const bool estimate_pi{true};
    static arma::dvec pi() { return arma::dvec{0.95, 0.05}; }
    static std::vector<std::string>
    prior_str(double nu, double s2, const arma::dvec& pi)
    {
        std::vector<std::string> prior_strings;
        prior_strings.reserve(4);
        prior_strings.emplace_back("BayesBπ");
        prior_strings.emplace_back("      ├─ αᵢ ~ (1-π) N(0, σ²ᵢ) + πδ₀");
        prior_strings.emplace_back(
            fmt::format("      ├─ {}", sigma_prior("ᵢ", nu, s2)));
        prior_strings.emplace_back(fmt::format("      └─ π = {}", pi.at(0)));
        return prior_strings;
    }
    static void sample(
        const bayes::GeneticEffect& design,
        bayes::GeneticEffectState& state,
        double* y_adj,
        arma::uvec& snp_tracker,
        double sigma_e,
        std::mt19937_64& rng)
    {
        GeneticTrait<BayesAlphabet::B>::sample(
            design, state, y_adj, snp_tracker, sigma_e, rng);
        state.pi.prop = dirichlet(state.pi.count + 1, rng);  // 1 for prior
    }
};

template <>
struct GeneticTrait<BayesAlphabet::C>
{
    static arma::dvec sigma(const arma::dmat& X) { return arma::dvec{0}; }

    static const bool estimate_pi{false};
    static arma::dvec pi() { return arma::dvec{0.95, 0.05}; }

    static std::vector<std::string>
    prior_str(double nu, double s2, const arma::dvec& pi)
    {
        std::vector<std::string> prior_strings;
        prior_strings.reserve(3);
        prior_strings.emplace_back("BayesC");
        prior_strings.emplace_back("      ├─ αᵢ ~ 0.05 N(0, σ²) + 0.95 δ₀");
        prior_strings.emplace_back(
            fmt::format("      └─ {}", sigma_prior("", nu, s2)));
        return prior_strings;
    }
    static void sample(
        const bayes::GeneticEffect& design,
        bayes::GeneticEffectState& state,
        double* y_adj,
        arma::uvec& snp_tracker,
        double sigma_e,
        std::mt19937_64& rng)
    {
        arma::dvec logpi = arma::log(state.pi.prop);

        arma::dvec& coeff = state.coeff;
        double* u = state.u.memptr();
        const arma::dmat& design_matrix = design.design_matrix;

        const double sigma = arma::as_scalar(state.sigma);

        double var_a{};
        std::normal_distribution<double> normal{0, 1};
        std::uniform_real_distribution<double> uniform{0, 1};
        const auto norm = static_cast<double>(design_matrix.n_rows) - 1;
        const int n = static_cast<int>(design_matrix.n_rows);

        for (size_t i = 0; i < coeff.n_elem; ++i)
        {
            const double old_i = coeff.at(i);
            const double* col_i = design_matrix.colptr(i);
            const double inv_scaler = 1 / (norm + sigma_e / sigma);

            double rhs = ddot_ptr(n, col_i, y_adj);
            if (old_i != 0)
            {
                rhs += norm * old_i;
            }
            double logdetV = log((sigma * norm / sigma_e) + 1);
            double uhat = rhs * inv_scaler;

            const double l_diff = (-0.5 * (logdetV - uhat * rhs / sigma_e))
                                  + logpi.at(1) - logpi.at(0);
            const double accept_prob = 1 / (1 + exp(l_diff));
            const bool has_effect = (uniform(rng) >= accept_prob);
            snp_tracker.at(i) = has_effect ? 1 : 0;

            double new_i{};
            if (has_effect)
            {
                new_i = normal(rng) * sqrt(sigma_e * inv_scaler) + uhat;
                const double diff = old_i - new_i;
                daxpy_ptr(n, diff, col_i, y_adj);
                daxpy_ptr(n, -diff, col_i, u);
                var_a += new_i * new_i;
            }
            else if (old_i != 0.0)
            {
                daxpy_ptr(n, old_i, col_i, y_adj);
                daxpy_ptr(n, -old_i, col_i, u);
            }
            coeff.at(i) = new_i;
        }
        state.pi.count.at(1) = arma::sum(snp_tracker);
        state.pi.count.at(0) = coeff.n_elem - state.pi.count.at(1);

        ScaledInvChiSq chi_squared{design.prior};
        chi_squared.update(var_a, static_cast<int>(state.pi.count.at(1)));
    }
};

template <>
struct GeneticTrait<BayesAlphabet::Cpi>
{
    static arma::dvec sigma(const arma::dmat& X)
    {
        return arma::zeros<arma::dvec>(1);
    }

    static const bool estimate_pi{true};
    static arma::dvec pi() { return arma::dvec{0.95, 0.05}; }
    static std::vector<std::string>
    prior_str(double nu, double s2, const arma::dvec& pi)
    {
        std::vector<std::string> prior_strings;
        prior_strings.reserve(4);
        prior_strings.emplace_back("BayesCπ");
        prior_strings.emplace_back("      ├─ αᵢ ~ (1-π) N(0, σ²) + πδ₀");
        prior_strings.emplace_back(
            fmt::format("      ├─ {}", sigma_prior("", nu, s2)));
        prior_strings.emplace_back(fmt::format("      └─ π = {}", pi.at(0)));
        return prior_strings;
    }
    static void sample(
        const bayes::GeneticEffect& design,
        bayes::GeneticEffectState& state,
        double* y_adj,
        arma::uvec& snp_tracker,
        double sigma_e,
        std::mt19937_64& rng)
    {
        GeneticTrait<BayesAlphabet::C>::sample(
            design, state, y_adj, snp_tracker, sigma_e, rng);
        state.pi.prop = dirichlet(state.pi.count + 1, rng);
    }
};

template <>
struct GeneticTrait<BayesAlphabet::R>
{
    static arma::dvec sigma(const arma::dmat& X)
    {
        return arma::dvec{0, 1e-4, 1e-3, 1e-2, 0};
    }  // this sigma from 0-3 is the scaler the last one is the ture sigma
    static const bool estimate_pi{true};
    static arma::dvec pi() { return arma::dvec{0.95, 0.02, 0.02, 0.01}; }
    static std::vector<std::string>
    prior_str(double nu, double s2, const arma::dvec& pi)
    {
        std::vector<std::string> prior_strings;
        prior_strings.reserve(2);
        prior_strings.emplace_back("BayesR");
        std::string pi_str = "      └─ αj ∼";
        std::string last_str = " (1";
        for (size_t i = 0; i < pi.n_elem; ++i)
        {
            if (i > 0)
            {
                pi_str += " +";
            }
            pi_str += fmt::format(" π{}N(0, {}σ²α)", i + 1, pi.at(i));
            last_str += fmt::format("-π{}", i + 1);
        }
        prior_strings.emplace_back(pi_str + last_str + ")δ₀");
        return prior_strings;
    }
    static void sample(
        const bayes::GeneticEffect& design,
        bayes::GeneticEffectState& state,
        double* y_adj,
        arma::uvec& snp_tracker,
        double sigma_e,
        std::mt19937_64& rng)
    {
        arma::dvec logpi = arma::log(state.pi.prop);
        arma::dvec& coeff = state.coeff;
        double* u = state.u.memptr();
        const arma::dmat& design_matrix = design.design_matrix;
        const auto norm = static_cast<double>(design_matrix.n_rows) - 1;

        for (size_t i = 0; i < state.coeff.n_elem; ++i)
        {
            // TODO: Implement R sampler
        }
    }
};

constexpr std::size_t to_index(BayesAlphabet e)
{
    return static_cast<std::size_t>(e);
}

using Fn = arma::dvec (*)(const arma::dmat&);
using FnPrior_str
    = std::vector<std::string> (*)(double, double, const arma::dvec&);
using FnSample = void (*)(
    const bayes::GeneticEffect&,
    bayes::GeneticEffectState&,
    double*,
    arma::uvec&,
    double,
    std::mt19937_64&);
using FnPi = arma::dvec (*)();

using BayesATrait = GeneticTrait<BayesAlphabet::A>;
using BayesRRTrait = GeneticTrait<BayesAlphabet::RR>;
using BayesBTrait = GeneticTrait<BayesAlphabet::B>;
using BayesBpiTrait = GeneticTrait<BayesAlphabet::Bpi>;
using BayesCTrait = GeneticTrait<BayesAlphabet::C>;
using BayesCpiTrait = GeneticTrait<BayesAlphabet::Cpi>;
using BayesRTrait = GeneticTrait<BayesAlphabet::R>;

inline constexpr std::array<Fn, to_index(BayesAlphabet::Count)>
    bayes_trait_sigma = {
        &BayesATrait::sigma,
        &BayesRRTrait::sigma,
        &BayesBTrait::sigma,
        &BayesBpiTrait::sigma,
        &BayesCTrait::sigma,
        &BayesCpiTrait::sigma,
        &BayesRTrait::sigma,
};

inline constexpr std::array<bool, to_index(BayesAlphabet::Count)>
    bayes_trait_estimate_pi = {
        BayesATrait::estimate_pi,
        BayesRRTrait::estimate_pi,
        BayesBTrait::estimate_pi,
        BayesBpiTrait::estimate_pi,
        BayesCTrait::estimate_pi,
        BayesCpiTrait::estimate_pi,
        BayesRTrait::estimate_pi,
};

inline constexpr std::array<FnSample, to_index(BayesAlphabet::Count)>
    bayes_trait_sample = {
        &BayesATrait::sample,
        &BayesRRTrait::sample,
        &BayesBTrait::sample,
        &BayesBpiTrait::sample,
        &BayesCTrait::sample,
        &BayesCpiTrait::sample,
        &BayesRTrait::sample,
};
inline constexpr std::array<FnPi, to_index(BayesAlphabet::Count)> bayes_trait_pi
    = {
        &BayesATrait::pi,
        &BayesRRTrait::pi,
        &BayesBTrait::pi,
        &BayesBpiTrait::pi,
        &BayesCTrait::pi,
        &BayesCpiTrait::pi,
        &BayesRTrait::pi,
};

inline constexpr std::array<FnPrior_str, to_index(BayesAlphabet::Count)>
    bayes_trait_prior_str = {
        &BayesATrait::prior_str,
        &BayesRRTrait::prior_str,
        &BayesBTrait::prior_str,
        &BayesBpiTrait::prior_str,
        &BayesCTrait::prior_str,
        &BayesCpiTrait::prior_str,
        &BayesRTrait::prior_str,
};
}  // namespace gelex
