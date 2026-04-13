/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef GELEX_ALGO_INFER_POSTERIOR_CALCULATOR_H_
#define GELEX_ALGO_INFER_POSTERIOR_CALCULATOR_H_

#include <Eigen/Core>

#include "gelex/types/posterior_summary.h"

namespace gelex::detail
{

class EigenThreadGuard
{
   public:
    EigenThreadGuard() : old_thread_count_(Eigen::nbThreads())
    {
        Eigen::setNbThreads(1);
    }

    ~EigenThreadGuard() { Eigen::setNbThreads(old_thread_count_); }

    EigenThreadGuard(const EigenThreadGuard&) = delete;
    EigenThreadGuard& operator=(const EigenThreadGuard&) = delete;
    EigenThreadGuard(EigenThreadGuard&&) = delete;
    EigenThreadGuard& operator=(EigenThreadGuard&&) = delete;

   private:
    int old_thread_count_;
};

namespace PosteriorCalculator
{

void compute_pve(
    PosteriorSummary& summary,
    const Eigen::Ref<const Eigen::VectorXd>& mean_coeffs,
    double phenotype_var);

}  // namespace PosteriorCalculator

}  // namespace gelex::detail

#endif  // GELEX_ALGO_INFER_POSTERIOR_CALCULATOR_H_
