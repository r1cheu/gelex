#pragma once
#include <cstddef>
#include <string>
#include <vector>

#include <armadillo>

#include "gelex/dist.h"
#include "gelex/model/effects.h"

namespace gelex
{

/**
 * @class Pi
 * @brief Structure to hold the proportion and counts of snps, use in BayesBpi,
 * BayesCpi, BayesR etc.
 */
struct Pi
{
    arma::dvec prop;
    arma::uvec count;
};

namespace bayes
{

/**
 * @class FixedEffect
 * @brief Structure to hold fixed effect information for Bayesian Models.
 *
 */
struct FixedEffect
{
    /**
     * @brief Initialize a FixedEffect with names, levels, and design matrix.
     *
     * @param names names of the fixed effects
     * @param levels levels for each fixed effect
     * @param design_matrix the design matrix for the fixed effects
     */
    FixedEffect(
        std::vector<std::string>&& names,
        std::vector<std::string>&& levels,
        arma::dmat&& design_matrix);

    arma::dmat design_matrix;
    arma::dvec cols_norm;
    std::vector<std::string> names;
    std::vector<std::string> levels;
};

struct FixedEffectState
{
    explicit FixedEffectState(size_t n_coeff);
    arma::dvec coeff;
};

struct RandomEffect
{
    RandomEffect(std::string&& name, arma::dmat&& design_matrix);

    std::string name;

    arma::dmat design_matrix;
    arma::dvec cols_norm;

    ScaledInvChiSqParams prior;
    arma::dvec sigma{0};  // set to dvec (not double) for consistency
};

struct RandomEffectState
{
    RandomEffectState(size_t n_coeff, const arma::dvec& init_sigma);
    arma::dvec coeff;
    arma::dvec sigma{0};  // set to dvec (not double) for consistency
};

struct GeneticEffect
{
    GeneticEffect(
        std::string&& name,
        arma::dmat&& design_matrix,
        BayesAlphabet type,
        arma::dvec&& sigma,
        arma::dvec&& pi);

    std::string name;
    arma::dmat design_matrix;

    ScaledInvChiSqParams prior;
    arma::dvec sigma;

    BayesAlphabet type;
    arma::dvec pi;

    arma::dvec mean;
    arma::dvec stddev;
};

struct GeneticEffectState
{
    GeneticEffectState(
        BayesAlphabet type_,
        size_t n_individual,
        size_t n_coeff,
        const arma::dvec& pi_prop,
        const arma::dvec& sigma_);

    BayesAlphabet type;

    arma::dvec coeff;
    arma::dvec u;

    Pi pi;
    double genetic_var{};
    double heritability{};
    arma::dvec sigma;
};

struct Residual
{
    std::string name{"e"};
    ScaledInvChiSqParams prior;
    double value{0.0};
};

using RandomEffectManager = Effects<RandomEffect>;
using GeneticEffectManager = Effects<GeneticEffect>;

std::vector<RandomEffectState> create_thread_states(
    const RandomEffectManager& designs);
std::vector<GeneticEffectState> create_thread_states(
    const GeneticEffectManager& designs);
}  // namespace bayes

template <typename Mat>
arma::dvec sum_square(const Mat& mat)
{
    arma::dvec result(mat.n_cols);

#pragma omp parallel for default(none) shared(result, mat)
    for (size_t i = 0; i < mat.n_cols; ++i)
    {
        result.at(i) = arma::dot(mat.col(i), mat.col(i));
    }
    return result;
}

arma::dvec compute_cols_var(const arma::dmat& mat);

}  // namespace gelex
