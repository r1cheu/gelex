#include <armadillo>

namespace chenx
{
using dmat = arma::mat;
using dvec = arma::vec;
void solver_chol(dmat& V, dmat& X, dvec& y);

}  // namespace chenx
