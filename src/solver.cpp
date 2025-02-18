#include "chenx/solver.h"
#include <armadillo>

namespace chenx
{
void solver_chol(dmat& V, dmat& X, dvec& y)
{
    char uplo = 'L';
    int n = static_cast<int>(V.n_cols);
    int info{};
    arma::lapack::potrf(&uplo, &n, V.memptr(), &n, &info);
    if (info != 0)
    {
        throw std::runtime_error("V Matrix is not symmetric positive definite");
    }

    dmat rhs = arma::join_horiz(X, y);
    int m = static_cast<int>(rhs.n_cols);
    arma::lapack::potrs(&uplo, &n, &m, V.memptr(), &n, rhs.memptr(), &n, &info);

    rhs.brief_print();
}
}  // namespace chenx
