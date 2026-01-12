#include "gelex/optim/constrain.h"

#include <Eigen/Core>

#include "gelex/logger.h"

namespace gelex
{

void constrain(Eigen::Ref<Eigen::VectorXd> varcmp, double y_variance)
{
    auto logger = logging::get();
    constexpr double constr_scale = 1e-6;
    const double limit = y_variance * constr_scale;

    Eigen::ArrayX<bool> mask = varcmp.array() < 0;
    auto num_constrained = mask.count();
    auto num_varcmp = varcmp.size();

    if (num_constrained == 0)
    {
        return;
    }
    if (num_constrained == varcmp.size())
    {
        logger->warn(
            "All variance components are constrained! The estimate is not "
            "reliable.");
        varcmp.fill(limit);
        return;
    }

    if (num_constrained > num_varcmp / 2)
    {
        logger->warn(
            "Half of the variance components are constrained! The estimate is "
            "not reliable.");
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
}
}  // namespace gelex
