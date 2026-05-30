#include "gelex/bayes/genetic/parameters.h"

#include <array>
#include <cmath>
#include <limits>
#include <optional>
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

MixtureProportion::MixtureProportion(Eigen::VectorXd initial_value)
    : initial_value_(std::move(initial_value))
{
    if (initial_value_.size() < 2)
    {
        throw GelexException(
            "MixtureProportion: initial_value must have at least 2 entries");
    }
    for (Eigen::Index i = 0; i < initial_value_.size(); ++i)
    {
        if (!std::isfinite(initial_value_(i)) || !(initial_value_(i) > 0))
        {
            throw GelexException(
                "MixtureProportion: initial_value entries must be finite and "
                "> 0");
        }
    }
    const double sum = initial_value_.sum();
    const double tol = std::numeric_limits<double>::epsilon()
                       * static_cast<double>(initial_value_.size()) * 64.0;
    if (std::abs(sum - 1.0) > tol)
    {
        throw GelexException(
            "MixtureProportion: initial_value entries must sum to 1");
    }
}

MixtureProportion::MixtureProportion(SimplexParameter parameter)
    : initial_value_(parameter.initial_value()),
      prior_(std::make_optional<DirichletPrior>(parameter.prior()))
{
}

}  // namespace gelex::bayes
