#pragma once
#include <string>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "../src/model/effects_manager.h"
#include "distribution.h"
#include "gelex/data/genotype_mmap.h"
#include "gelex/model/effects.h"

namespace gelex
{
namespace bayes
{
/**
 * @class Pi
 * @brief Structure to hold the proportion and counts of snps, use in BayesBpi,
 * BayesCpi, BayesR etc.
 */
struct Pi
{
    Eigen::VectorXd prop;
    Eigen::VectorXi count;
};

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
        Eigen::MatrixXd&& design_matrix);

    Eigen::MatrixXd design_matrix;
    Eigen::VectorXd cols_norm;

    std::vector<std::string> names;
};

struct FixedStatus
{
    explicit FixedStatus(Eigen::Index n_coeff);
    Eigen::VectorXd coeff;
};

struct RandomEffect
{
    RandomEffect(std::string&& name, Eigen::MatrixXd&& design_matrix);

    std::string name;
    std::vector<std::string> levels;

    Eigen::MatrixXd design_matrix;
    Eigen::VectorXd cols_norm;

    detail::ScaledInvChiSqParams prior;

    Eigen::VectorXd effect_variance{
        {0}};  // set to dvec (not double) for consistency
};

struct RandomStatus
{
    explicit RandomStatus(const RandomEffect& effect);

    Eigen::VectorXd coeff;
    Eigen::VectorXd effect_variance{
        {0}};  // set to dvec (not double) for consistency
};

struct AdditiveEffect
{
    AdditiveEffect(
        std::string&& name,
        GenotypeMap&& design_matrix,
        Eigen::VectorXd&& sigma,
        Eigen::VectorXd&& pi);

    std::string name;
    GenotypeMap design_matrix;
    Eigen::VectorXd cols_norm;

    detail::ScaledInvChiSqParams prior;
    Eigen::VectorXd marker_variance;

    Eigen::VectorXd pi;

    bool is_monomorphic(Eigen::Index snp_index) const;
    Eigen::Index num_mono() const;
};

struct AdditiveStatus
{
    explicit AdditiveStatus(const AdditiveEffect& effect);

    Eigen::VectorXd coeff;
    Eigen::VectorXd u;
    Eigen::VectorXi tracker;

    Pi pi;
    double effect_variance{};
    double heritability{};
    Eigen::VectorXd marker_variance;
};

struct DominantEffect
{
    DominantEffect(
        std::string&& name,
        GenotypeMap&& design_matrix,
        double prior_mean,
        double prior_var);

    std::string name;
    GenotypeMap design_matrix;
    Eigen::VectorXd cols_norm;

    double prior_mean{};
    double prior_var{};

    bool is_monomorphic(Eigen::Index snp_index) const;
    Eigen::Index num_mono() const;
};

struct DominantStatus
{
    explicit DominantStatus(const DominantEffect& effect);

    Eigen::VectorXd coeff;
    Eigen::VectorXd u;

    double effect_variance{};
    double heritability{};
};

struct Residual
{
    std::string name{"e"};
    detail::ScaledInvChiSqParams prior;
    double value{0.0};
};

using RandomEffectManager = detail::Effects<RandomEffect>;

std::vector<RandomStatus> create_chain_states(
    const RandomEffectManager& effects);

}  // namespace bayes
}  // namespace gelex
