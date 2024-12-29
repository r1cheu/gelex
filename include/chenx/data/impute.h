#pragma once
#include <armadillo>

namespace chenx
{
using namespace arma;

dvec MeanImpute(dmat& genotype);
dvec MedianImpute(dmat& genotype);
void ValueImpute(dmat& genotype, const dvec& values);
}  // namespace chenx
