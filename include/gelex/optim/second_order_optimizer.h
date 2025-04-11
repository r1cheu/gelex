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
    dvec step(const GBLUP& model) override;

   private:
    virtual dvec computeFirstGrad(const GBLUP& model);

    virtual double compute_hess_elem(const GBLUP& model, uword i, uword j) = 0;
    dmat compute_hess(const GBLUP& model);
};

class NewtonRaphsonOptimizer : public SecondOrderOptimizer
{
    using SecondOrderOptimizer::SecondOrderOptimizer;

   public:
    std::string name() const noexcept override { return "NewtonRaphson"; }

   private:
    double compute_hess_elem(const GBLUP& model, uword i, uword j) override;
};

class FisherScoringOptimizer : public SecondOrderOptimizer
{
    using SecondOrderOptimizer::SecondOrderOptimizer;

   public:
    std::string name() const noexcept override { return "FisherScoring"; }

   private:
    double compute_hess_elem(const GBLUP& model, uword i, uword j) override;
};

class AverageInformationOptimizer : public SecondOrderOptimizer
{
    using SecondOrderOptimizer::SecondOrderOptimizer;

   public:
    std::string name() const noexcept override { return "AverageInformation"; }

   private:
    double compute_hess_elem(const GBLUP& model, uword i, uword j) override;
};

}  // namespace gelex
