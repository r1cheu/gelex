#pragma once
#include <cstddef>
#include <memory>
#include <string>

#include <spdlog/logger.h>
#include <armadillo>

#include <chenx/utils.h>
#include "chenx/model/linear_mixed_model.h"

namespace chenx
{
class OptimizerBase
{
   public:
    explicit OptimizerBase(size_t max_iter = 20, double tol = 1e-8)
        : max_iter_{max_iter}, tol_{tol}, logger_{Logger::logger()}
    {
    }
    OptimizerBase(const OptimizerBase&) = default;
    OptimizerBase(OptimizerBase&&) = delete;
    OptimizerBase& operator=(const OptimizerBase&) = default;
    OptimizerBase& operator=(OptimizerBase&&) = delete;
    virtual ~OptimizerBase() = default;

    virtual std::string name() const noexcept = 0;

    size_t max_iter() const noexcept { return max_iter_; }
    double tol() const noexcept { return tol_; }
    bool converged() const noexcept { return converged_; }
    double obj_func_diff() const noexcept { return obj_func_diff_; }

    void set_max_iter(size_t max_iter) noexcept { max_iter_ = max_iter; }
    void set_tol(double tol) noexcept { tol_ = tol; }
    void set_converged(bool converged) noexcept { converged_ = converged; }

    virtual bool Optimize(LinearMixedModel& model);
    virtual dvec Step(const LinearMixedModel& model) = 0;

   protected:
    dvec Constrain(dvec sigma, double y_var);

   private:
    size_t max_iter_;
    double tol_;
    std::shared_ptr<spdlog::logger> logger_;

    bool converged_{};
    dvec old_param_;
    double old_obj_func_value_{};
    double obj_func_diff_{};

    double VecDiff(const dvec& new_param);
    double ObjFunctionDiff(double new_value);
    void CheckConvergence(const dvec& new_param, double new_value);
};
}  // namespace chenx
