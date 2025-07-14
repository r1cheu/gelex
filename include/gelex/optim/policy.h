#pragma once
#include "gelex/optim/optimizer.h"

namespace gelex
{
struct ExpectationMaximizationPolicy
{
    static arma::dvec apply(Optimizer& optim, GBLUP& model);
};

struct AverageInformationPolicy
{
    static arma::dvec apply(Optimizer& optim, GBLUP& model);
};

}  // namespace gelex
