// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_REML_REML_BUFFER_H_
#define GELEX_FREQ_REML_REML_BUFFER_H_

#include <Eigen/Core>

#include "gelex/freq/model.h"

namespace gelex
{

class RemlBuffer
{
   public:
    explicit RemlBuffer(const gelex::FreqModel& model);

    auto phenotype_variance() const -> double { return phenotype_variance_; }
    auto num_individuals() const -> Eigen::Index { return num_individuals_; }

    // tr(P) = tr(V^{-1}) - tr(XtViX_inv * ViX' * ViX)
    [[nodiscard]] auto trace_proj() const -> double;

    // tr(P * K) = tr(V^{-1} * K) - tr(XtViX_inv * ViX' * K * ViX)
    [[nodiscard]] auto trace_proj_k(
        const Eigen::Ref<const Eigen::MatrixXd>& K) const -> double;

    // computed matrices, public for Policy access
    // Projection is represented lazily: P = V^{-1} - ViX * XtViX_inv * ViX'
    Eigen::MatrixXd V;
    Eigen::VectorXd Py;
    Eigen::MatrixXd ViX;
    Eigen::MatrixXd XtViX_inv;
    double logdet_v{};
    double logdet_xvx{};

    // for AI policy
    Eigen::MatrixXd hess_inv;
    Eigen::MatrixXd dvpy;        // n x n_comp, each column is K_i * Py
    Eigen::VectorXd first_grad;  // first derivative

   private:
    Eigen::Index num_individuals_{};
    double phenotype_variance_{};
};

}  // namespace gelex

#endif  // GELEX_FREQ_REML_REML_BUFFER_H_
