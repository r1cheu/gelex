#pragma once

#include <armadillo>

namespace chenx
{
using arma::dmat;

dmat AdditiveGrm(dmat& genotype);
void AddChunkGrm(dmat& genotype, dmat& grm);
dmat DomainanceGrm(dmat& genotype);
dmat ComputeGRM(dmat& genotype);

}  // namespace chenx
