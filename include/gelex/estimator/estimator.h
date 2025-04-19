#pragma once
#include <cstddef>
#include <memory>
#include <string_view>

#include <armadillo>

#include <fmt/ranges.h>
#include <spdlog/logger.h>
#include "gelex/model/gblup.h"
#include "gelex/optim/base_optimizer.h"

namespace gelex
{

class Estimator
{
   public:
    Estimator(std::string_view optimizer, size_t max_iter, double tol);
    void set_optimizer(std::string_view optimizer, double tol);
    void fit(GBLUP& model, bool em_init = true, bool verbose = true);

   private:
    dvec compute_beta(GBLUP& model);
    dmat compute_u(GBLUP& model);
    std::unique_ptr<OptimizerBase> optimizer_;
    std::shared_ptr<spdlog::logger> logger_;

    uint64_t max_iter_{};
    double tol_{};
    std::string optimizer_name_;
    bool converged_{};
};

inline dvec compute_h2(const GBLUP& model)
{
    return model.genetic_sigma() / arma::sum(model.sigma());
}

template <typename eT>
auto blue_vec(const arma::Col<eT>& vec)
{
    return fmt::styled(fmt::join(vec, ", "), fmt::fg(fmt::color::blue_violet));
}

};  // namespace gelex
