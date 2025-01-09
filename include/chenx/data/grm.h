#pragma once

#include <armadillo>

namespace chenx
{
using arma::dmat;
void HandleNaN(dmat& genotype);
void Normalize(dmat& genotype);
dmat AddGrm(dmat& genotype);
void AddGrmChunk(dmat& genotype, dmat& grm);
dmat DomGrm(dmat& genotype);
void DomGrmChunk(dmat& genotype, dmat& grm);
dmat ComputeGRM(dmat& genotype);

}  // namespace chenx
