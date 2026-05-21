#include "gelex/model/bayes/prior_specs.h"

#include <cmath>

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

VarianceSpec::VarianceSpec(
    PositiveScalar<double> initial_value,
    ScaledInvChiSqPrior prior)
    : initial_value_(initial_value), prior_(prior)
{
}

VarianceSpec::VarianceSpec(double initial_value, ScaledInvChiSqPrior prior)
    : VarianceSpec(PositiveScalar<double>{initial_value}, prior)
{
}

}  // namespace gelex::bayes
