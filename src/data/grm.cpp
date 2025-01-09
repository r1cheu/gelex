#include "chenx/data/grm.h"

#include <armadillo>
#include <cmath>

namespace chenx
{
using arma::dvec;

void HandleNaN(dmat& genotype)
{
    if (genotype.has_nan())
    {
        genotype.replace(arma::datum::nan, 1.0);
    }
}

void Normalize(dmat& genotype)
{
    dvec pA = mean(genotype, 0) / 2;
    dvec pa = 1 - pA;
    genotype.replace(2.0, 0.0);
    genotype.each_row() -= 2 * (pA % pa);
}

dmat ComputeGRM(dmat& genotype)
{
    dmat grm = genotype * genotype.t();
    return grm / trace(grm) * static_cast<double>(grm.n_rows);
}

dmat AddGrm(dmat& genotype)
{
    HandleNaN(genotype);
    genotype.each_row() -= mean(genotype, 0);
    return ComputeGRM(genotype);
}

void AddGrmChunk(dmat& genotype, dmat& grm)
{
    HandleNaN(genotype);
    genotype.each_row() -= mean(genotype, 0);
    grm += genotype * genotype.t();
}

dmat DomGrm(dmat& genotype)
{
    HandleNaN(genotype);
    Normalize(genotype);
    return ComputeGRM(genotype);
}

void DomGrmChunk(dmat& genotype, dmat& grm)
{
    HandleNaN(genotype);
    Normalize(genotype);
    grm += genotype * genotype.t();
}

}  // namespace chenx
