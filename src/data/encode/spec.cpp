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

#include "gelex/data/encode/spec.h"

#include "gelex/data/encode/types.h"
#include "gelex/data/genotype_method.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

auto encoding_spec_from_method(GeneticMode effect, GenotypeMethod method)
    -> EncodingSpec
{
    EncodingSpec spec;
    spec.effect = effect;
    spec.normalization = is_center(method) ? Normalization::Center
                                           : Normalization::CenterScale;
    spec.moment_basis
        = is_hwe(method) ? MomentBasis::Theoretical : MomentBasis::Empirical;

    if (is_noia(method))
    {
        spec.dominance_code = DominanceCode::NOIA;
    }
    else if (is_orthogonal(method))
    {
        spec.dominance_code = DominanceCode::HWE;
    }
    else
    {
        spec.dominance_code = DominanceCode::Het;
    }

    return spec;
}

}  // namespace gelex
