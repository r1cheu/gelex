#include "chenx/data/grm.h"

#include "armadillo"

namespace chenx
{
dmat ComputeGRM(dmat& genotype)
{
    dmat grm = genotype.t() * genotype;
    return grm / trace(grm) * static_cast<double>(grm.n_cols);
}

dmat AdditiveGrm(dmat& genotype)
{
    genotype.each_col() -= mean(genotype, 1);
    return ComputeGRM(genotype);
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
