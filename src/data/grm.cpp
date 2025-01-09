#include "chenx/data/grm.h"

#include <armadillo>
#include <cmath>

namespace chenx
{
using arma::dvec;
dmat ComputeGRM(dmat& genotype)
{
    if (genotype.has_nan())
    {
        genotype.replace(arma::datum::nan, 1.0);
    }
    dmat grm = genotype * genotype.t();
    return grm / trace(grm) * static_cast<double>(grm.n_rows);
}

dmat AdditiveGrm(dmat& genotype)
{
    if (genotype.has_nan())
    {
        genotype.replace(arma::datum::nan, 1.0);
    }

    genotype.each_row() -= mean(genotype, 0);
    return ComputeGRM(genotype);
}

void AddChunkGrm(dmat& genotype, dmat& grm)
{
    if (genotype.has_nan())
    {
        genotype.replace(arma::datum::nan, 1.0);
    }

    genotype.each_row() -= mean(genotype, 0);
    grm += genotype * genotype.t();
}

dmat DomainanceGrm(dmat& genotype)
{
    dvec pA = mean(genotype, 1) / 2;
    dvec pa = 1 - pA;
    genotype.replace(2.0, 0.0);
    genotype.each_col() -= 2 * (pA % pa);
    return ComputeGRM(genotype);
}

}  // namespace chenx
