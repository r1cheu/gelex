#pragma once
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string_view>

#include <armadillo>

#include <fmt/ranges.h>
#include <spdlog/logger.h>
#include "gelex/model/effects.h"
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
    void report_variance(
        const std::string& category,
        const std::vector<size_t>& indices,
        const GBLUP& model);
    uint64_t max_iter_{};
    double tol_{};
    std::string optimizer_name_;
    bool converged_{};
};

template <typename eT>
auto blue_vec(const arma::Col<eT>& vec)
{
    return fmt::styled(fmt::join(vec, ", "), fmt::fg(fmt::color::blue_violet));
}

std::string join_formula(
    const std::vector<uint64_t>& indices,
    const GroupEffectManager& effects,
    std::string_view sep);

std::string join_name(
    const std::vector<uint64_t>& indices,
    const GroupEffectManager& effects,
    std::string_view sep);

};  // namespace gelex
