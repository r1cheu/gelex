#include "chenx/model/linear_mixed_model.h"

#include <stdexcept>
#include <vector>

#include <armadillo>

namespace chenx
{
LinearMixedModel::LinearMixedModel(
    dvec&& y,
    dmat&& X,
    dcube&& covar_matrices_rand,
    std::vector<std::string>&& rand_names)
    : y_{std::move(y)},
      X_{std::move(X)},
      zkzt_{covar_matrices_rand},
      rand_names_{std::move(rand_names)}
{
    y_var_ = var(y_);
    uword n{X_.n_rows}, m{X_.n_cols};
    // zkztr_ = ComputeZKZtR(std::move(covar_matrices_rand));
    uword n_rands = zkzt_.n_slices + 1;
    pdv_.set_size(n, n, n_rands);
    set_sigma(dvec(
        n_rands, arma::fill::value(y_var_ / static_cast<double>(n_rands))));
    rand_names_.emplace_back("Residual");
}

double LinearMixedModel::ComputeLogLikelihood() const
{
    return -0.5
           * (logdet_v_ + log_det_sympd(tx_vinv_x_)
              + as_scalar(y_.t() * proj_y_));
}

/*dmat LinearMixedModel::ComputeZKZ(const sp_dmat& z, const dmat& k)*/
/*{*/
/*    bool z_identity = CheckIdentity(z);*/
/*    bool k_identity = CheckIdentity(k);*/
/**/
/*    if (k_identity)*/
/*    {*/
/*        return z_identity ? k : dmat(z * z.t());*/
/*    }*/
/*    return z_identity ? k : dmat(z * k * z.t());*/
/*}*/
/**/
/*dcube LinearMixedModel::ComputeZKZtR(*/
/*    dcube&& covar_matrices_rand)  // the matrices will be moved in and
 * release*/
/*                                  // the memory*/
/*{*/
/*    dcube rands{std::move(covar_matrices_rand)};*/
/*    auto n = rands.n_rows;*/
/*    auto n_rands = rands.n_slices;*/
/*    dcube result(n, n, n_rands + 1, arma::fill::zeros);*/
/**/
/*    for (size_t i = 0; i < n_rands; ++i)*/
/*    {*/
/*        result.slice(i) = rands.slice(i);*/
/*    }*/
/*    result.slice(n_rands).eye();*/
/**/
/*    return result;*/
/*}*/
/**/
void LinearMixedModel::ComputeV()
{
    v_ = sigma_.at(0) * zkzt_.slice(0);
    for (size_t i{1}; i < zkzt_.n_slices; ++i)
    {
        v_ += sigma_.at(i) * zkzt_.slice(i);
    }

    v_.diag() += sigma_.back();
}

void LinearMixedModel::ComputeProj()
{
    logdet_v_ = VinvLogdet(v_);  // from here the v_ matrix is inverted
    dmat vinv_x = v_ * X_;
    tx_vinv_x_ = X_.t() * vinv_x;

    proj_
        = v_
          - vinv_x
                * solve(tx_vinv_x_, vinv_x.t(), arma::solve_opts::likely_sympd);
    proj_y_ = proj_ * y_;
}

void LinearMixedModel::ComputePdV()
{
    for (size_t i = 0; i < zkzt_.n_slices; ++i)
    {
        pdv_.slice(i) = proj_ * zkzt_.slice(i);
    }

    pdv_.slice(zkzt_.n_slices) = proj_;
}

double LinearMixedModel::VinvLogdet(dmat& V)
{
    char uplo = 'L';
    int n = static_cast<int>(V.n_cols);
    int info{};
    arma::lapack::potrf(&uplo, &n, V.memptr(), &n, &info);
    if (info != 0)
    {
        throw std::runtime_error("V Matrix is not symmetric positive definite");
    }
    double logdet = accu(2.0 * log(diagvec(V)));
    arma::lapack::potri(&uplo, &n, V.memptr(), &n, &info);

    if (info != 0)
    {
        throw std::runtime_error("Error during inverse V matrix");
    }
    V = arma::symmatl(V);
    return logdet;
}

}  // namespace chenx
