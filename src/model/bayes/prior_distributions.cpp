#include "gelex/model/bayes/prior_distributions.h"

#include <cmath>
#include <utility>

#include <Eigen/Core>

#include "gelex/exception.h"

namespace gelex::bayes
{

ScaledInvChiSqPrior::ScaledInvChiSqPrior(
    double degrees_of_freedom,
    double scale)
    : degrees_of_freedom_(degrees_of_freedom), scale_(scale)
{
    const bool is_proper = degrees_of_freedom_ > 0 && scale_ > 0
                           && std::isfinite(degrees_of_freedom_)
                           && std::isfinite(scale_);
    const bool is_flat
        = degrees_of_freedom_ == -2 && scale_ == 0;  // flat prior

    if (!is_proper && !is_flat)
    {
        throw GelexException(
            "ScaledInvChiSqPrior: prior must be either proper "
            "(df > 0, scale > 0, both finite) "
            "or flat sentinel (df = -2, scale = 0)");
    }
}

DirichletPrior::DirichletPrior(Eigen::VectorXd concentration)
    : concentration_(std::move(concentration))
{
    if (concentration_.size() < 2)
    {
        throw GelexException(
            "DirichletPrior: concentration must have at least 2 entries");
    }
    for (Eigen::Index i = 0; i < concentration_.size(); ++i)
    {
        if (!std::isfinite(concentration_(i)) || !(concentration_(i) > 0))
        {
            throw GelexException(
                "DirichletPrior: concentration entries must be finite and "
                "> 0");
        }
    }
}

}  // namespace gelex::bayes
