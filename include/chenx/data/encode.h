#pragma once
#include <armadillo>

namespace chenx
{
using namespace arma;
dmat ComputeHybirdValue(const dmat& genotype, const dmat& phenotype);
void HybridEncode(dmat& genotype, const dmat& hybird_value);
}  // namespace chenx
