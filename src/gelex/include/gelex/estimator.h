#pragma once
#include <cstddef>
#include <memory>
#include <string_view>

#include <armadillo>

#include <spdlog/logger.h>
#include "gelex/model/linear_mixed_model.h"
#include "gelex/optim/base_optimizer.h"

namespace gelex
{

class Estimator
{
   public:
    Estimator(std::string_view optimizer, size_t max_iter, double tol);
    void set_optimizer(std::string_view optimizer, size_t max_iter, double tol);
    void Fit(LinearMixedModel& model, bool em_init = true, bool verbose = true);

   private:
    template <typename OptimizerType>
    static std::unique_ptr<OptimizerBase> CreateOptimizer(
        size_t max_iter,
        double tol)
    {
        return std::make_unique<OptimizerType>(max_iter, tol);
    }
    static dvec ComputeBeta(LinearMixedModel& model);
    static dmat ComputeU(LinearMixedModel& model);
    std::unique_ptr<OptimizerBase> optimizer_;
    std::shared_ptr<spdlog::logger> logger_;
};
};  // namespace gelex
