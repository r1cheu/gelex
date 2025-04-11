#pragma once
#include <armadillo>
#include <string>

#include "gelex/model/gblup.h"
#include "gelex/optim/base_optimizer.h"

namespace gelex
{
class ExpectationMaximizationOptimizer : public OptimizerBase
{
    using OptimizerBase::OptimizerBase;

   public:
    dvec Step(const GBLUP& model) override;

    std::string name() const noexcept override
    {
        return "ExpectationMaximization";
    }
};
}  // namespace gelex
