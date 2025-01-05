#pragma once
#include <armadillo>
#include <string>

#include "chenx/model/linear_mixed_model.h"
#include "chenx/optim/base_optimizer.h"

namespace chenx
{
class ExpectationMaximizationOptimizer : public OptimizerBase
{
    using OptimizerBase::OptimizerBase;

   public:
    dvec Step(const LinearMixedModel& model) override;

    std::string name() const noexcept override
    {
        return "ExpectationMaximization";
    }
};
}  // namespace chenx
