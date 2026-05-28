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

#include "gelex/algo/infer/fixed_samples.h"

#include "gelex/model/bayes/designs.h"
#include "gelex/model/bayes/state.h"
#include "gelex/types/fixed_designs.h"

namespace gelex
{

FixedSamples::FixedSamples(const FixedDesign& design)
    : names(design.names), levels(design.levels), n_coeffs_(design.X.cols())
{
}

RandomSamples::RandomSamples(const bayes::RandomDesign& design)
    : levels(design.levels), n_coeffs_(design.X.cols())
{
}

void FixedSamples::store(const bayes::FixedState& state)
{
    coeffs_stats_.update(state.coeffs);
}

void RandomSamples::store(const bayes::RandomState& state)
{
    coeffs_stats_.update(state.coeffs);
    variance_stats_.update(state.variance);
}

}  // namespace gelex
