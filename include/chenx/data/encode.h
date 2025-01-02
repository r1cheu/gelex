#pragma once
#include <armadillo>

namespace chenx
{
using arma::dmat;
using arma::dvec;
dmat ComputeHybirdValue(const dmat& genotype, const dvec& phenotype);
void HybridEncode(dmat& genotype, const dmat& hybird_value);
}  // namespace chenx
