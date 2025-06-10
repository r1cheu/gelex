#pragma once

#include <armadillo>
#include <utility>

#include "gelex/model/gblup.h"
#include "gelex/optim/base_optimizer.h"

namespace gelex
{
class ExpectationMaximizationOptimizer : public OptimizerBase
{
    using OptimizerBase::OptimizerBase;

   public:
    explicit ExpectationMaximizationOptimizer(OptimizerBase&& base)
        : OptimizerBase(std::move(base))
    {
    }

   private:
    void step_inner(GBLUP& model) override;
};

class AverageInformationOptimizer : public OptimizerBase
{
    using OptimizerBase::OptimizerBase;

   public:
    explicit AverageInformationOptimizer(OptimizerBase&& base)
        : OptimizerBase(std::move(base))
    {
    }

   private:
    void step_inner(GBLUP& model) override;
    arma::dmat compute_hess();
};

}  // namespace gelex
