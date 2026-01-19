#ifndef GELEX_GWAS_ASSOCIATION_TEST_H
#define GELEX_GWAS_ASSOCIATION_TEST_H

#include <Eigen/Core>
#include "gelex/types/assoc_input.h"

namespace gelex::gwas
{

// Perform Wald test for a single SNP
// residual: y - Xβ̂ (phenotype minus fixed effects)
// v_inv: pre-computed V⁻¹ from null model
void wald_test(AssocInput& input, AssocOutput& output);
}  // namespace gelex::gwas

#endif  // GELEX_GWAS_ASSOCIATION_TEST_H
