#pragma once
#include "chenx/optim/base_optimizer.h"

namespace chenx
{
class SecondOrderOptimizer : public OptimizerBase
{
    using OptimizerBase::OptimizerBase;

   public:
    dvec Step(const LinearMixedModel& model);

   private:
    virtual dvec ComputeFirstGrad(const LinearMixedModel& model);

    virtual double
    ComputeHessElement(const LinearMixedModel& model, uword i, uword j)
        = 0;
    dmat ComputeHess(const LinearMixedModel& model);
};

class NewtonRaphsonOptimizer : public SecondOrderOptimizer
{
    using SecondOrderOptimizer::SecondOrderOptimizer;

   private:
    double ComputeHessElement(const LinearMixedModel& model, uword i, uword j)
        override;
};

class FisherScoringOptimizer : public SecondOrderOptimizer
{
    using SecondOrderOptimizer::SecondOrderOptimizer;

   private:
    double ComputeHessElement(const LinearMixedModel& model, uword i, uword j)
        override;
};

class AverageInformationOptimizer : public SecondOrderOptimizer
{
    using SecondOrderOptimizer::SecondOrderOptimizer;

   private:
    double ComputeHessElement(const LinearMixedModel& model, uword i, uword j)
        override;
};

}  // namespace chenx
