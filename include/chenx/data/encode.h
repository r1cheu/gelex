#pragma once
#include <armadillo>

namespace chenx
{
using arma::dmat;
dmat ComputeHybirdValue(const dmat& genotype, const dmat& phenotype);
void HybridEncode(dmat& genotype, const dmat& hybird_value);
}  // namespace chenx
