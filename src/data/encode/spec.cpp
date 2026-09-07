// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/data/encode/spec.h"

#include "gelex/data/encode/types.h"
#include "gelex/data/genotype_method.h"
#include "gelex/genetic_mode.h"

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
