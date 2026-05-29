#include "gelex/model/bayes/prior_parameters.h"

#include <array>
#include <cmath>
#include <limits>
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

MarkerVariance::MarkerVariance(
    MarkerVarianceLayout layout,
    VarianceParameter parameter)
    : layout_(layout), parameter_(parameter)
{
}

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

MixtureProportion::MixtureProportion(
    SimplexParameter parameter,
    UpdatePolicy update)
    : parameter_(std::move(parameter)), update_(update)
{
}

}  // namespace gelex::bayes
