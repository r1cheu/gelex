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

#ifndef GELEX_ESTIMATOR_BAYES_PARAMS_H_
#define GELEX_ESTIMATOR_BAYES_PARAMS_H_
#include <cstddef>
#include <stdexcept>

#include <Eigen/Core>

namespace gelex
{
struct MCMCParams
{
    MCMCParams(Eigen::Index n_iters, Eigen::Index n_burnin, Eigen::Index n_thin)
        : n_iters{n_iters},
          n_burnin{n_burnin},
          n_thin{n_thin},
          n_records{(n_iters - n_burnin) / n_thin}
    {
        if (n_burnin >= n_iters)
        {
            throw std::invalid_argument(
                "n_burnin must be smaller than n_iters");
        }
    }
    Eigen::Index n_iters;
    Eigen::Index n_burnin;
    Eigen::Index n_thin;
    Eigen::Index n_records;
};
}  // namespace gelex

#endif  // GELEX_ESTIMATOR_BAYES_PARAMS_H_
