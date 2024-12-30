#pragma once
#include <armadillo>

namespace chenx
{
using arma::dmat;
using arma::dvec;
dvec MeanImpute(dmat& genotype);
dvec MedianImpute(dmat& genotype);
void ValueImpute(dmat& genotype, const dvec& values);
}  // namespace chenx
