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

#ifndef GELEX_ESTIMATOR_FREQ_REML_H_
#define GELEX_ESTIMATOR_FREQ_REML_H_

#include <memory>

#include "gelex/data/data_pipe.h"

namespace gelex
{
struct AssocInput;
class SampleManager;

auto load_data_for_reml(const DataPipe::Config& config) -> DataPipe;

auto reml(
    const DataPipe::Config& config,
    size_t max_iter = 100,
    double tol = 1e-8,
    bool em_init = true,
    bool verbose = true) -> std::
    tuple<std::shared_ptr<SampleManager>, Eigen::MatrixXd, Eigen::VectorXd>;

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_FREQ_REML_H_
