// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/freq/reml/reml_buffer.h"

#include <Eigen/Core>

#include "gelex/freq/model.h"

namespace gelex
{

RemlBuffer::RemlBuffer(const FreqModel& model)
    : num_individuals_(model.num_individuals()),
      phenotype_variance_(model.phenotype_variance())
{
    const auto n_fixed = model.fixed().X().cols();
    V.resize(num_individuals_, num_individuals_);
    Py.resize(num_individuals_);
    ViX.resize(num_individuals_, n_fixed);
    XtViX_inv.resize(n_fixed, n_fixed);

    // preallocate for AI policy
    // n_comp = 1 (residual) + n_random
    auto n_comp = static_cast<Eigen::Index>(1 + model.random().size());
    dvpy.resize(num_individuals_, n_comp);
    first_grad.resize(n_comp);
}

auto RemlBuffer::trace_proj() const -> double
{
    return V.trace() - XtViX_inv.cwiseProduct(ViX.transpose() * ViX).sum();
}

auto RemlBuffer::trace_proj_k(const Eigen::Ref<const Eigen::MatrixXd>& K) const
    -> double
{
    return V.cwiseProduct(K).sum()
           - XtViX_inv.cwiseProduct(ViX.transpose() * K * ViX).sum();
}

}  // namespace gelex
