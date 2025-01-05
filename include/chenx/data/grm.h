#pragma once

#include <armadillo>

namespace chenx
{
using arma::dmat;

dmat AdditiveGrm(dmat& genotype);
dmat DomainanceGrm(dmat& genotype);
dmat ComputeGRM(dmat& genotype);

}  // namespace chenx
