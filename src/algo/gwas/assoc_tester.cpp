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

#include "gelex/algo/gwas/assoc_tester.h"

#include "gelex/algo/gwas/joint_tester.h"
#include "gelex/algo/gwas/single_tester.h"
#include "gelex/exception.h"

namespace gelex
{

auto AssocTester::make(
    AssocTestType type,
    GeneticMode mode,
    GenotypeProcessMethod geno_method) -> std::unique_ptr<AssocTester>
{
    if (!geno_method.is_center())
    {
        throw GelexException(
            "assoc --geno-method supports only center-family methods: "
            "CH(center-hwe), OCH(orth-center-hwe), C(center), "
            "OC(orth-center)");
    }

    switch (type)
    {
        case AssocTestType::Single:
            return std::make_unique<SingleTester>(mode, geno_method);
        case AssocTestType::Joint:
            return std::make_unique<JointTester>(geno_method);
    }
    throw GelexException("Unknown AssocTestType");
}

}  // namespace gelex
