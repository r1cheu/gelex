#pragma once

namespace gelex
{

// Prior distribution constants for Bayesian models
namespace prior_constants
{

// Inverse-Gamma prior parameters for marker variance
// Shape parameter for marker variance prior
constexpr double MARKER_VARIANCE_SHAPE = 4.0;
// Scale parameter multiplier for marker variance prior
constexpr double MARKER_VARIANCE_SCALE_MULTIPLIER = 0.5;

// Inverse-Gamma prior parameters for residual variance
// Shape parameter for residual variance prior
constexpr double RESIDUAL_VARIANCE_SHAPE = -2.0;
// Scale parameter for residual variance prior
constexpr double RESIDUAL_VARIANCE_SCALE = 0.0;

// Inverse-Gamma prior parameters for random effects variance
// Shape parameter for random effects variance prior
constexpr double RANDOM_EFFECTS_SHAPE = 4.0;
// Scale parameter for random effects variance prior
constexpr double RANDOM_EFFECTS_SCALE = 0.0;

// Default non-zero marker proportion for non mixture model
constexpr double NON_MIXTURE_PROPORTION = 1.0;
}  // namespace prior_constants

}  // namespace gelex
