#include <armadillo>

#include "gelex/dist.h"
#include "gelex/estimator/gibbs/base.h"
#include "gelex/model/bayes_effects.h"

namespace gelex
{

template <BayesAlphabet>
struct GeneticTrait;

static constexpr double DEFAULT_SIGMA = 0.01;
template <>
struct GeneticTrait<BayesAlphabet::A>
{
    static arma::dvec sigma(const arma::dmat& X)
    {
        return arma::dvec(X.n_cols, arma::fill::value(DEFAULT_SIGMA));
    }
    static arma::dvec pi() { return arma::dvec(); }
    static void sample(
        const GeneticEffectDesign& design,
        GeneticEffectState& state,
        double* y_adj,
        arma::uvec& snp_tracker,
        double sigma_e,
        std::mt19937_64& rng)
    {
        arma::dvec& coeff = state.coeff;
        double* u = state.u.memptr();

        const arma::dvec& cols_norm = design.cols_norm;
        const arma::dmat& design_mat = design.design_mat;
        arma::dvec& sigma = state.sigma;
        const arma::dvec& cols_var = design.cols_var;

        std::normal_distribution<double> normal{0, 1};
        const int n = static_cast<int>(design_mat.n_rows);

        for (size_t i = 0; i < coeff.n_elem; ++i)
        {
            if (cols_var.at(i) == 0)
            {
                continue;
            }
            const double old_i = coeff.at(i);
            const double new_sigma = sample_scale_inv_chi_squared(
                rng,
                design.prior.nu + 1,
                ((old_i * old_i) + (design.prior.s2 * design.prior.nu)));
            const double* col_i = design_mat.colptr(i);
            const double norm = cols_norm.at(i);

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
    static arma::dvec sigma(const arma::dmat& X)
    {
        return arma::dvec{DEFAULT_SIGMA};
    }
    static arma::dvec pi() { return arma::dvec(); }
    static void sample(
        const GeneticEffectDesign& design,
        GeneticEffectState& state,
        double* y_adj,
        arma::uvec& snp_tracker,
        double sigma_e,
        std::mt19937_64& rng)
    {
        arma::dvec& coeff = state.coeff;
        double* u = state.u.memptr();
        const arma::dvec& cols_norm = design.cols_norm;
        const arma::dmat& design_mat = design.design_mat;
        const double sigma_g = arma::as_scalar(state.sigma);
        const arma::dvec& cols_var = design.cols_var;

        const double sigma_e_sqrt = sqrt(sigma_e);
        const double inv_scaler_base = sigma_e / sigma_g;

        std::normal_distribution<double> normal{0, 1};
        const int n = static_cast<int>(design_mat.n_rows);

        for (size_t idx = 0; idx < coeff.n_elem; ++idx)
        {
            if (cols_var.at(idx) == 0)
            {
                continue;
            }

            const double old_i = coeff.at(idx);
            const double norm = cols_norm.at(idx);
            const double inv_scaler = 1.0 / (norm + inv_scaler_base);

            const double* col_i = design_mat.colptr(idx);
            double rhs = ddot_ptr(n, col_i, y_adj) + (norm * old_i);
            const double new_i = (normal(rng) * sigma_e_sqrt * sqrt(inv_scaler))
                                 + (rhs * inv_scaler);

            coeff.at(idx) = new_i;
            const double diff = old_i - new_i;
            daxpy_ptr(n, diff, col_i, y_adj);
            daxpy_ptr(n, -diff, col_i, u);
        }
        const double ssq = arma::dot(state.coeff, state.coeff);
        state.sigma.at(0) = sample_scale_inv_chi_squared(
            rng,
            design.prior.nu
                + static_cast<double>(coeff.n_elem - design.n_zero_var_snp),
            (ssq + (design.prior.s2 * design.prior.nu)));
    }
};

template <>
struct GeneticTrait<BayesAlphabet::B>
{
    static arma::dvec sigma(const arma::dmat& X)
    {
        return arma::zeros<arma::dvec>(X.n_cols);
    }

    static arma::dvec pi() { return arma::dvec{0.95, 0}; }
    static void sample(
        const GeneticEffectDesign& design,
        GeneticEffectState& state,
        double* y_adj,
        arma::uvec& snp_tracker,
        double sigma_e,
        std::mt19937_64& rng)
    {
        arma::dvec logpi = arma::log(state.pi.prop);

        arma::dvec& coeff = state.coeff;
        double* u = state.u.memptr();
        const arma::dvec& cols_norm = design.cols_norm;
        const arma::dmat& design_mat = design.design_mat;
        arma::dvec& sigma = state.sigma;
        const arma::dvec& cols_var = design.cols_var;

        std::normal_distribution<double> normal{0, 1};
        std::uniform_real_distribution<double> uniform{0, 1};

        const int n = static_cast<int>(design_mat.n_rows);

        for (size_t i = 0; i < coeff.n_elem; ++i)
        {
            if (cols_var.at(i) == 0)
            {
                continue;
            }

            const double old_i = coeff.at(i);
            const double new_sigma = sample_scale_inv_chi_squared(
                rng,
                design.prior.nu + 1,
                ((old_i * old_i) + (design.prior.s2 * design.prior.nu)));

            const double* col_i = design_mat.colptr(i);
            const double norm = cols_norm.at(i);
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
        state.pi.count.at(0)
            = coeff.n_elem - design.n_zero_var_snp - state.pi.count.at(1);
    }
};

template <>
struct GeneticTrait<BayesAlphabet::Bpi>
{
    static arma::dvec sigma(const arma::dmat& X)
    {
        return arma::zeros<arma::dvec>(X.n_cols);
    }
    static arma::dvec pi() { return arma::dvec{0.95, 0}; }
    static void sample(
        const GeneticEffectDesign& design,
        GeneticEffectState& state,
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
    static arma::dvec sigma(const arma::dmat& X)
    {
        return arma::dvec{DEFAULT_SIGMA};
    }
    static arma::dvec pi() { return arma::dvec{0.95, 0.05}; }
    static void sample(
        const GeneticEffectDesign& design,
        GeneticEffectState& state,
        double* y_adj,
        arma::uvec& snp_tracker,
        double sigma_e,
        std::mt19937_64& rng)
    {
        arma::dvec logpi = arma::log(state.pi.prop);

        arma::dvec& coeff = state.coeff;
        double* u = state.u.memptr();
        const arma::dvec& cols_norm = design.cols_norm;
        const arma::dmat& design_mat = design.design_mat;
        const arma::dvec& cols_var = design.cols_var;

        const double sigma = arma::as_scalar(state.sigma);

        double var_a{};
        std::normal_distribution<double> normal{0, 1};
        std::uniform_real_distribution<double> uniform{0, 1};

        const int n = static_cast<int>(design_mat.n_rows);
        for (size_t i = 0; i < coeff.n_elem; ++i)
        {
            if (cols_var.at(i) == 0)
            {
                continue;
            }
            const double old_i = coeff.at(i);
            const double* col_i = design_mat.colptr(i);
            const double norm = cols_norm.at(i);
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
        state.pi.count.at(0)
            = coeff.n_elem - design.n_zero_var_snp - state.pi.count.at(1);

        state.sigma.at(0) = sample_scale_inv_chi_squared(
            rng,
            design.prior.nu + static_cast<double>(state.pi.count.at(1)),
            (var_a + (design.prior.s2 * design.prior.nu)));
    }
};

template <>
struct GeneticTrait<BayesAlphabet::Cpi>
{
    static arma::dvec sigma(const arma::dmat& X)
    {
        return arma::zeros<arma::dvec>(1);
    }
    static arma::dvec pi() { return arma::dvec{0.95, 0}; }
    static void sample(
        const GeneticEffectDesign& design,
        GeneticEffectState& state,
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
    static arma::dvec pi() { return arma::dvec{0.95, 0.02, 0.02, 0.01}; }
    static void sample(
        const GeneticEffectDesign& effect,
        GeneticEffectState& state,
        double* y_adj,
        arma::uvec& snp_tracker,
        double sigma_e,
        std::mt19937_64& rng)
    {
        for (size_t i = 0; i < state.coeff.n_elem; ++i)
        {
            if (effect.cols_var.at(i) == 0)
            {
                continue;
            }
            // TODO: Implement R sampler
        }
    }
};

constexpr std::size_t to_index(BayesAlphabet e)
{
    return static_cast<std::size_t>(e);
}

using Fn = arma::dvec (*)(const arma::dmat&);
using FnSample = void (*)(
    const GeneticEffectDesign&,
    GeneticEffectState&,
    double*,
    arma::uvec&,
    double,
    std::mt19937_64&);
using FnPi = arma::dvec (*)();

inline constexpr std::array<Fn, to_index(BayesAlphabet::Count)>
    bayes_trait_sigma = {
        &GeneticTrait<BayesAlphabet::A>::sigma,
        &GeneticTrait<BayesAlphabet::RR>::sigma,
        &GeneticTrait<BayesAlphabet::B>::sigma,
        &GeneticTrait<BayesAlphabet::Bpi>::sigma,
        &GeneticTrait<BayesAlphabet::C>::sigma,
        &GeneticTrait<BayesAlphabet::Cpi>::sigma,
        &GeneticTrait<BayesAlphabet::R>::sigma,
};

inline constexpr std::array<FnSample, to_index(BayesAlphabet::Count)>
    bayes_trait_sample = {
        &GeneticTrait<BayesAlphabet::A>::sample,
        &GeneticTrait<BayesAlphabet::RR>::sample,
        &GeneticTrait<BayesAlphabet::B>::sample,
        &GeneticTrait<BayesAlphabet::Bpi>::sample,
        &GeneticTrait<BayesAlphabet::C>::sample,
        &GeneticTrait<BayesAlphabet::Cpi>::sample,
        &GeneticTrait<BayesAlphabet::R>::sample,
};
inline constexpr std::array<FnPi, to_index(BayesAlphabet::Count)> bayes_trait_pi
    = {
        &GeneticTrait<BayesAlphabet::A>::pi,
        &GeneticTrait<BayesAlphabet::RR>::pi,
        &GeneticTrait<BayesAlphabet::B>::pi,
        &GeneticTrait<BayesAlphabet::Bpi>::pi,
        &GeneticTrait<BayesAlphabet::C>::pi,
        &GeneticTrait<BayesAlphabet::Cpi>::pi,
        &GeneticTrait<BayesAlphabet::R>::pi,
};
};  // namespace gelex
