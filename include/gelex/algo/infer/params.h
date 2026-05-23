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

#ifndef GELEX_ALGO_INFER_PARAMS_H_
#define GELEX_ALGO_INFER_PARAMS_H_

#include <Eigen/Core>

namespace gelex
{

namespace mcmc
{

struct Params
{
    Eigen::Index n_iters{};
    Eigen::Index n_burn_in{};
    Eigen::Index n_thin{};
    Eigen::Index checkpoint_step{};

    auto n_records() const -> Eigen::Index
    {
        return (n_iters - n_burn_in) / n_thin;
    }
};

}  // namespace mcmc

}  // namespace gelex

#endif  // GELEX_ALGO_INFER_PARAMS_H_
