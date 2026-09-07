// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/freq/gwas/assoc_tester.h"

#include <memory>

#include "gelex/data/genotype_method.h"
#include "gelex/exception.h"
#include "gelex/freq/gwas/assoc_type.h"
#include "gelex/freq/gwas/joint_tester.h"
#include "gelex/freq/gwas/single_tester.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

auto AssocTester::make(
    AssocType type,
    GeneticMode mode,
    GenotypeMethod geno_method) -> std::unique_ptr<AssocTester>
{
    if (!is_center(geno_method))
    {
        throw GelexException(
            "assoc --geno-method supports only center-family methods: "
            "CH(center-hwe), OCH(orth-center-hwe), C(center), "
            "OC(orth-center)");
    }

    switch (type)
    {
        case AssocType::Single:
            return std::make_unique<SingleTester>(mode, geno_method);
        case AssocType::Joint:
            return std::make_unique<JointTester>(geno_method);
    }
    throw GelexException("Unknown AssocType");
}

}  // namespace gelex
