#include "gelex/bayes/genetic/parameters.h"

#include <array>
#include <cmath>
#include <limits>
#include <utility>

#include <Eigen/Core>

#include "gelex/exception.h"

namespace gelex::bayes
{

SharedMarkerVariance::SharedMarkerVariance(VarianceParameter parameter)
    : parameter_(parameter)
{
}

PerMarkerVariance::PerMarkerVariance(VarianceParameter parameter)
    : parameter_(parameter)
{
}

JointSharedMarkerVariance::JointSharedMarkerVariance(
    std::array<SharedMarkerVariance, 2> variances)
    : variances_(std::move(variances))
{
}

FixedProportion::FixedProportion(Eigen::VectorXd value)
    : value_(std::move(value))
{
    if (value_.size() < 2)
    {
        throw GelexException(
            "FixedMixtureProportion: value must have at least 2 entries");
    }
    for (Eigen::Index i = 0; i < value_.size(); ++i)
    {
        if (!std::isfinite(value_(i)) || !(value_(i) > 0))
        {
            throw GelexException(
                "FixedMixtureProportion: value entries must be finite and "
                "> 0");
        }
    }
    const double sum = value_.sum();
    const double tol = std::numeric_limits<double>::epsilon()
                       * static_cast<double>(value_.size()) * 64.0;
    if (std::abs(sum - 1.0) > tol)
    {
        throw GelexException(
            "FixedMixtureProportion: value entries must sum to 1");
    }
}

SampledProportion::SampledProportion(SimplexParameter parameter)
    : parameter_(std::move(parameter))
{
}

}  // namespace gelex::bayes
