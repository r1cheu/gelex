#ifndef GELEX_TYPES_FREQ_EFFECT_H_
#define GELEX_TYPES_FREQ_EFFECT_H_

#include <string>
#include <vector>

#include <Eigen/Core>

#include "../src/types/fixed_effects.h"

namespace gelex::freq
{

enum class GrmType : uint8_t
{
    A,
    D,
    AD,
    AA,
    DD,
    Unknown
};

struct RandomEffect
{
    std::string name;
    std::vector<std::string> levels;
    Eigen::MatrixXd K;  // covariance matrix (n x n)
};

struct GeneticEffect
{
    GrmType type;
    Eigen::MatrixXd K;  // GRM matrix
};

struct FixedState
{
    explicit FixedState(const gelex::FixedEffect& effect);
    Eigen::VectorXd coeff;
    Eigen::VectorXd se;
};

struct RandomState
{
    explicit RandomState(const RandomEffect& effect);
    std::string name;
    Eigen::VectorXd blup;
    double variance{};
    double variance_se{};
};

struct GeneticState
{
    explicit GeneticState(const GeneticEffect& effect);
    GrmType type;
    Eigen::VectorXd ebv;
    double variance{};
    double variance_se{};
    double heritability{};
    double heritability_se{};
};

struct ResidualState
{
    double variance{};
    double variance_se{};
};

}  // namespace gelex::freq
#endif  // GELEX_TYPES_FREQ_EFFECT_H_
