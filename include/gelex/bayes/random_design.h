// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_RANDOM_DESIGN_H_
#define GELEX_BAYES_RANDOM_DESIGN_H_

#include <Eigen/Core>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/data/dataframe/key_type.h"

namespace gelex
{

template <KeyType Key>
class DataFrame;

}  // namespace gelex

namespace gelex::bayes
{

class RandomDesign
{
   public:
    RandomDesign(
        std::string name,
        std::vector<std::string> column_names,
        Eigen::MatrixXd X);

    RandomDesign(const RandomDesign&) = default;
    auto operator=(const RandomDesign&) -> RandomDesign& = default;
    RandomDesign(RandomDesign&&) noexcept = default;
    auto operator=(RandomDesign&&) noexcept -> RandomDesign& = default;
    ~RandomDesign() = default;

    [[nodiscard]] auto name() const noexcept -> std::string_view
    {
        return name_;
    }

    [[nodiscard]] auto column_names() const noexcept
        -> std::span<const std::string>
    {
        return column_names_;
    }

    [[nodiscard]] auto X() const noexcept -> const Eigen::MatrixXd&
    {
        return matrix_;
    }

    [[nodiscard]] auto xtx_diag() const noexcept -> const Eigen::VectorXd&
    {
        return xtx_diag_;
    }

   private:
    std::string name_;
    std::vector<std::string> column_names_;
    Eigen::MatrixXd matrix_;
    Eigen::VectorXd xtx_diag_;
};

[[nodiscard]] auto make_random_designs(const DataFrame<std::string>& frame)
    -> std::vector<RandomDesign>;

[[nodiscard]] auto make_quantitative_random_design(
    const DataFrame<std::string>& frame,
    std::string name) -> RandomDesign;

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_RANDOM_DESIGN_H_
