#pragma once
#include <string>

#include <armadillo>

#include "gelex/model/gblup.h"
#include "gelex/optim/base_optimizer.h"

namespace gelex
{
class SecondOrderOptimizer : public OptimizerBase
{
    using OptimizerBase::OptimizerBase;

   public:
    dvec Step(const GBLUP& model) override;

   private:
    virtual dvec ComputeFirstGrad(const GBLUP& model);

    virtual double ComputeHessElement(const GBLUP& model, uword i, uword j) = 0;
    dmat ComputeHess(const GBLUP& model);
};

class NewtonRaphsonOptimizer : public SecondOrderOptimizer
{
    using SecondOrderOptimizer::SecondOrderOptimizer;

   public:
    std::string name() const noexcept override { return "NewtonRaphson"; }

   private:
    double ComputeHessElement(const GBLUP& model, uword i, uword j) override;
};

class FisherScoringOptimizer : public SecondOrderOptimizer
{
    using SecondOrderOptimizer::SecondOrderOptimizer;

   public:
    std::string name() const noexcept override { return "FisherScoring"; }

   private:
    double ComputeHessElement(const GBLUP& model, uword i, uword j) override;
};

class AverageInformationOptimizer : public SecondOrderOptimizer
{
    using SecondOrderOptimizer::SecondOrderOptimizer;

   public:
    std::string name() const noexcept override { return "AverageInformation"; }

   private:
    double ComputeHessElement(const GBLUP& model, uword i, uword j) override;
};

}  // namespace gelex
