// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_MARKER_COVARIATE_H
#define GELEX_BAYES_MARKER_COVARIATE_H

#include <Eigen/Core>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/dataframe/key_type.h"

namespace gelex
{
template <KeyType Key>
class DataFrame;
}  // namespace gelex

namespace gelex::bayes
{
class MarkerCovariate;

auto make_marker_covariate(
    DataFrame<std::string> frame,
    const DataFrame<std::string>& marker_metadata) -> MarkerCovariate;

class MarkerCovariate
{
   public:
    MarkerCovariate(const MarkerCovariate&) = delete;
    auto operator=(const MarkerCovariate&) -> MarkerCovariate& = delete;
    MarkerCovariate(MarkerCovariate&&) noexcept = default;
    auto operator=(MarkerCovariate&&) noexcept -> MarkerCovariate& = default;
    ~MarkerCovariate() = default;

    [[nodiscard]] auto annotation_names() const noexcept
        -> std::span<const std::string>
    {
        return annotation_names_;
    }

    [[nodiscard]] auto X() const noexcept -> const Eigen::MatrixXd&
    {
        return values_;
    }

   private:
    MarkerCovariate(
        std::vector<std::string> annotation_names,
        Eigen::MatrixXd values)
        : annotation_names_(std::move(annotation_names)),
          values_(std::move(values))
    {
    }

    std::vector<std::string> annotation_names_;
    Eigen::MatrixXd values_;

    friend auto make_marker_covariate(
        DataFrame<std::string> frame,
        const DataFrame<std::string>& marker_metadata) -> MarkerCovariate;
};
}  // namespace gelex::bayes

#endif  // GELEX_BAYES_MARKER_COVARIATE_H
