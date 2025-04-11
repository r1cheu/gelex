#include <fmt/format.h>
#include <armadillo>

#include "gelex/python/dense_caster.h"
#include "gelex/python/sparse_caster.h"

namespace bind
{
arma::sp_dmat to_sparse(
    iarr1d indice,
    iarr1d indptr,
    arr1d values,
    uint64_t row,
    uint64_t col)
{
    arma::uvec indptr_vec = arma::conv_to<arma::uvec>::from(to_arma(indptr));
    arma::uvec indice_vec = arma::conv_to<arma::uvec>::from(to_arma(indice));

    return arma::sp_dmat(
        arma::conv_to<arma::uvec>::from(indice_vec),
        arma::conv_to<arma::uvec>::from(indptr_vec),
        to_arma(values),
        row,
        col,
        false);
}

}  // namespace bind
