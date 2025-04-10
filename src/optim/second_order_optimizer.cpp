#include "gelex/optim/second_order_optimizer.h"

#include <armadillo>

namespace gelex
{
dvec SecondOrderOptimizer::Step(const GBLUP& model)
{
    dvec sigma{model.sigma()};
    uword num_random_effects{sigma.n_elem};
    dvec first_grad = ComputeFirstGrad(model);
    dmat hess = ComputeHess(model);
    dmat hess_inv;
    if (!pinv(hess_inv, hess))
    {
        throw std::runtime_error("Hessian matrix is not invertible!");
    }
    dvec delta = -hess_inv * first_grad;
    if (obj_func_diff() > 1)
    {
        sigma += 0.316 * delta;
    }
    else
    {
        sigma += delta;
    }
    return OptimizerBase::Constrain(sigma, model.y_var());
};

dvec SecondOrderOptimizer::ComputeFirstGrad(const GBLUP& model)
{
    const dvec& sigma = model.sigma();
    const dvec& y = model.y();
    const dvec& proj_y = model.proj_y();
    const dcube& pdv = model.pdv();
    uword n = sigma.n_elem;

    dvec first_grad(n, arma::fill::zeros);

    for (size_t i{0}; i < n; ++i)
    {
        first_grad.at(i) = -0.5
                           * (trace(pdv.slice(i))
                              - as_scalar(y.t() * pdv.slice(i) * proj_y));
    }
    return first_grad;
};

dmat SecondOrderOptimizer::ComputeHess(const GBLUP& model)
{
    uword n = model.sigma().n_elem;
    dmat hess(n, n, arma::fill::zeros);
    for (size_t i{0}; i < n; ++i)
    {
        for (size_t j{i}; j < n; ++j)
        {
            hess.at(i, j) = ComputeHessElement(model, i, j);
            if (i != j)
            {
                hess.at(j, i) = hess.at(i, j);
            }
        }
    }
    return hess;
}

double
NewtonRaphsonOptimizer::ComputeHessElement(const GBLUP& model, uword i, uword j)
{
    const dmat& pdvi = model.pdv().slice(i);
    const dmat& pdvj = model.pdv().slice(j);

    return 0.5 * trace(pdvi * pdvj)
           - as_scalar(model.y().t() * pdvi * pdvj * model.proj_y());
}

double
FisherScoringOptimizer::ComputeHessElement(const GBLUP& model, uword i, uword j)
{
    return -0.5 * trace(model.pdv().slice(i) * model.pdv().slice(j));
}

double AverageInformationOptimizer::ComputeHessElement(
    const GBLUP& model,
    uword i,
    uword j)
{
    return -0.5
           * as_scalar(
               model.y().t() * model.pdv().slice(i) * model.pdv().slice(j)
               * model.proj_y());
}
};  // namespace gelex
