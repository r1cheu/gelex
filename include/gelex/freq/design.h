// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_DESIGN_H_
#define GELEX_FREQ_DESIGN_H_

#include <Eigen/Core>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace gelex
{
class FixedDesign;
}

namespace gelex::freq
{

enum class RandomKind : std::uint8_t
{
    Discrete,      // one-hot ZZ^T from a discrete factor
    Quantitative,  // linear kernel ZZ^T from a quantitative matrix
    Grm,           // genomic relationship matrix; chromosome-partitionable
    Interaction,   // Hadamard product of two base kernels, rescaled
};

struct RandomDesign
{
    std::string name;
    std::optional<std::vector<std::string>> levels;
    std::optional<Eigen::MatrixXd> Z;  // skip if identity
    Eigen::MatrixXd K;                 // Kernels
    RandomKind kind = RandomKind::Discrete;
};

struct RandomState
{
    explicit RandomState(const RandomDesign& design);
    RandomState() = default;
    Eigen::VectorXd blup;  // sample-level random predictions
    double variance{};
    double variance_se{};
    double variance_ratio{};
    double variance_ratio_se{};
    bool at_boundary{};  // clamped to the constraint floor; Wald test invalid
};

struct FixedState
{
    explicit FixedState(const gelex::FixedDesign& design);
    Eigen::VectorXd coeffs;
    Eigen::VectorXd se;
};

struct ResidualState
{
    double variance{};
    double variance_se{};
};

}  // namespace gelex::freq

#endif  // GELEX_FREQ_DESIGN_H_
