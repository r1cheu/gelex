#pragma once
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../src/types/fixed_effects.h"

namespace gelex::freq
{

struct RandomEffect
{
    std::string name;
    std::vector<std::string> levels;
    Eigen::SparseMatrix<double> Z;  // incidence matrix
    Eigen::SparseMatrix<double> K;  // covariance matrix
};

struct GeneticEffect
{
    std::string name;
    Eigen::MatrixXd K;  // GRM matrix
};

struct FixedState
{
    explicit FixedState(const FixedEffect& effect);
    Eigen::VectorXd coeff;
    Eigen::VectorXd se;
};

struct RandomState
{
    explicit RandomState(const RandomEffect& effect);
    std::string name;
    Eigen::VectorXd blup;
    double variance;
};

struct GeneticState
{
    explicit GeneticState(const GeneticEffect& effect);
    std::string name;
    Eigen::VectorXd ebv;
    double variance;
};

struct ResidualState
{
    double variance{};
    double heritability{};
};

}  // namespace gelex::freq
