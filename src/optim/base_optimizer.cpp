#include "chenx/optim/base_optimizer.h"

#include "chenx/model/linear_mixed_model.h"

namespace chenx
{

bool OptimizerBase::Optimize(LinearMixedModel& model)
{
    for (size_t i{0}; i < max_iter(); ++i)
    {
        model.set_sigma(Step(model));
        double loglike{model.ComputeLogLikelihood()};
        CheckConvergence(model.sigma(), loglike);
        if (converged())
        {
            break;
            return true;
        }
    }
    return false;
}

double OptimizerBase::VecDiff(const dvec& new_param)
{
    if (old_param_.is_empty())
    {
        old_param_ = new_param;
        return 1e10;
    }
    double diff = norm(new_param - old_param_) / norm(new_param);
    old_param_ = new_param;
    return diff;
}

double OptimizerBase::ObjFunctionDiff(double new_value)
{
    if (old_obj_func_value_ == 0)
    {
        old_obj_func_value_ = new_value;
        return 1e10;
    }
    double diff{new_value - old_obj_func_value_};
    old_obj_func_value_ = new_value;
    return diff;
}

void OptimizerBase::CheckConvergence(const dvec& new_param, double new_value)
{
    double param_diff = VecDiff(new_param);
    double obj_func_value = ObjFunctionDiff(new_value);
    bool negative = obj_func_value < 0.0;
    obj_func_value = std::abs(obj_func_value);

    if (param_diff < tol_
        && (obj_func_value < 1e-4 || (negative && obj_func_value < 1e-2)))
    {
        converged_ = true;
    }
}

dvec OptimizerBase::Constrain(dvec sigma, double y_var)
{
    constexpr double constr_scale = 1e-6;
    arma::vec constrained_sigma = sigma;
    arma::uvec constrained = arma::find(sigma < 0.0);
    if (constrained.is_empty())
    {
        return constrained_sigma;
    }
    arma::uvec unconstrained = arma::find(sigma > 0.0);
    constrained_sigma(constrained).fill(y_var * constr_scale);

    double delta = arma::accu(y_var * constr_scale - sigma(constrained))
                   / static_cast<double>(unconstrained.n_elem);

    for (const auto& i : unconstrained)
    {
        if (sigma(i) > delta)
        {
            constrained_sigma.at(i) -= delta;
        }
    }

    if (constrained.n_elem > sigma.n_elem / 2)
    {
        std::cerr << "Half of the variance components are constrained! The "
                     "estimate is not reliable.\n";
    }

    return constrained_sigma;
}

}  // namespace chenx
