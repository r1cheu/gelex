#include <armadillo>

#include "dense_caster.h"

namespace bind
{
arma::sp_dmat to_sparse(
    iarr1d indice,
    iarr1d indptr,
    arr1d values,
    uint64_t row,
    uint64_t col);

}  // namespace bind
