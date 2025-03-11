#include "gelex/optim/base_optimizer.h"

#include <cstdlib>

#include <fmt/color.h>
#include <fmt/ranges.h>
#include <armadillo>

#include "gelex/model/linear_mixed_model.h"
#include "gelex/utils.h"

namespace gelex
{

bool OptimizerBase::Optimize(LinearMixedModel& model)
{
    double time_cost{};
    for (size_t i{1}; i <= max_iter(); ++i)
    {
        {
            Timer time(time_cost);
            model.set_sigma(Step(model));
        }
        double loglike{model.ComputeLogLikelihood()};
        CheckConvergence(model.sigma(), loglike);
        logger_->info(
            "Iter {:d}: logL={:.2f}, varcomp=[{:.4f}] {:.3f}s",
            i,
            loglike,
            fmt::styled(
                fmt::join(model.sigma(), " "),
                fmt::fg(fmt::color::blue_violet)),
            time_cost);
        if (converged())
        {
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
    double obj_func_diff = ObjFunctionDiff(new_value);
    bool negative = obj_func_diff < 0.0;
    obj_func_diff = std::abs(obj_func_diff);
    obj_func_diff_ = obj_func_diff;

    if (param_diff < tol_
        && (obj_func_diff < 1e-4 || (negative && obj_func_diff < 1e-2)))
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

}  // namespace gelex
