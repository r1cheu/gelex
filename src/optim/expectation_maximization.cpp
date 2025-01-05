#include "chenx/optim/expectation_maximization.h"

#include <armadillo>

#include "chenx/model/linear_mixed_model.h"
#include "chenx/optim/base_optimizer.h"

namespace chenx
{

dvec ExpectationMaximizationOptimizer::Step(const LinearMixedModel& model)
{
    dvec sigma{model.sigma()};
    dvec sigma2{sigma % sigma};
    const dvec& y = model.y();
    const dcube& pdv = model.pdv();
    uword n = y.n_elem;

    for (size_t i{0}; i < model.sigma().n_elem; ++i)
    {
        sigma.at(i) = as_scalar(
                          sigma2.at(i) * (y.t() * pdv.slice(i) * model.proj_y())
                          + trace(
                              sigma.at(i) * arma::eye(n, n)
                              - sigma2.at(i) * pdv.slice(i)))
                      / static_cast<double>(n);
    }
    return OptimizerBase::Constrain(sigma, model.y_var());
}

}  // namespace chenx
