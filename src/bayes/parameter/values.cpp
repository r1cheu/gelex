#include "gelex/bayes/parameter/values.h"

#include <cmath>
#include <limits>
#include <utility>

#include <Eigen/Core>

#include "gelex/exception.h"

namespace gelex::bayes
{

VarianceParameter::VarianceParameter(
    double initial_value,
    ScaledInvChiSqPrior prior)
    : initial_value_(initial_value), prior_(prior)
{
    if (!std::isfinite(initial_value_) || !(initial_value_ > 0))
    {
        throw GelexException(
            "VarianceParameter: initial_value must be finite and positive");
    }
}

SimplexParameter::SimplexParameter(
    Eigen::VectorXd initial_value,
    DirichletPrior prior)
    : initial_value_(std::move(initial_value)), prior_(std::move(prior))
{
    if (initial_value_.size() < 2)
    {
        throw GelexException(
            "SimplexParameter: initial_value must have at least 2 entries");
    }
    for (Eigen::Index i = 0; i < initial_value_.size(); ++i)
    {
        if (!std::isfinite(initial_value_(i)) || !(initial_value_(i) > 0))
        {
            throw GelexException(
                "SimplexParameter: initial_value entries must be finite and "
                "> 0");
        }
    }
    const double sum = initial_value_.sum();
    const double tol = std::numeric_limits<double>::epsilon()
                       * static_cast<double>(initial_value_.size()) * 64.0;
    if (std::abs(sum - 1.0) > tol)
    {
        throw GelexException(
            "SimplexParameter: initial_value entries must sum to 1");
    }
    if (prior_.size() != initial_value_.size())
    {
        throw GelexException(
            "SimplexParameter: initial_value and prior must have "
            "the same size");
    }
}

ProbabilityParameter::ProbabilityParameter(
    double initial_value,
    BetaPrior prior)
    : initial_value_(initial_value), prior_(prior)
{
    if (!std::isfinite(initial_value_) || !(initial_value_ > 0)
        || !(initial_value_ < 1))
    {
        throw GelexException(
            "ProbabilityParameter: initial_value must be finite and in "
            "(0, 1)");
    }
}

}  // namespace gelex::bayes
