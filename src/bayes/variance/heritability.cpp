// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/variance/heritability.h"

#include <Eigen/Core>
#include <fmt/format.h>

#include "gelex/exception.h"

namespace gelex
{

auto heritability_draws(
    const Eigen::Ref<const Eigen::RowVectorXd>& explained,
    const Eigen::Ref<const Eigen::RowVectorXd>& total_genetic,
    const Eigen::Ref<const Eigen::RowVectorXd>& residual) -> Eigen::RowVectorXd
{
    if (explained.size() != total_genetic.size()
        || explained.size() != residual.size())
    {
        throw GelexException(
            fmt::format(
                "heritability: draw counts differ (explained {}, total {}, "
                "residual {})",
                explained.size(),
                total_genetic.size(),
                residual.size()));
    }
    return explained.array() / (total_genetic.array() + residual.array());
}

}  // namespace gelex
