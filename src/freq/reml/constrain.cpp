// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/freq/reml/constrain.h"

#include <Eigen/Core>

namespace gelex
{

auto constrain(Eigen::Ref<Eigen::VectorXd> varcmp, double y_variance)
    -> Eigen::ArrayX<bool>
{
    const double limit = constraint_floor(y_variance);

    Eigen::ArrayX<bool> mask = varcmp.array() < 0;
    auto num_constrained = mask.count();
    auto num_varcmp = varcmp.size();

    if (num_constrained == 0)
    {
        return mask;
    }
    if (num_constrained == varcmp.size())
    {
        varcmp.fill(limit);
        return mask;
    }

    double delta = 0.0;
    for (int i = 0; i < num_varcmp; ++i)
    {
        if (mask[i])
        {
            delta += limit - varcmp[i];
            varcmp[i] = limit;
        }
    }

    delta /= static_cast<double>(num_varcmp - num_constrained);
    for (int i = 0; i < num_varcmp; ++i)
    {
        if (!mask[i] && varcmp[i] > delta)
        {
            varcmp[i] -= delta;
        }
    }
    return mask;
}
}  // namespace gelex
