#include "chenx/model/linear_mixed_model.h"
#include <spdlog/spdlog.h>

#include "chenx/logger.h"
#include "chenx/utils.h"

namespace chenx
{
LinearMixedModel::LinearMixedModel(
    dvec&& y,
    dmat&& X,
    std::vector<sp_dmat>&& Z,
    dcube&& covar_matrices_rand,
    std::vector<std::string>&& rand_names)
    : y_{std::move(y)},
      X_{std::move(X)},
      Z_{std::move(Z)},
      covar_matrices_rand_{std::move(covar_matrices_rand)},
      rand_names_{std::move(rand_names)},
      y_var_{var(y)},
      logger_{spdlog::get(logger_name)}
{
    uword n{X.n_rows}, m{X.n_cols};
    zkztr_ = ComputeZKZtR();
    uword n_rands = zkztr_.n_slices;
    beta_.zeros(m);
    sigma_.set_size(n_rands).fill(var(y) / static_cast<double>(n_rands));

    proj_y_.zeros(n);
    v_.zeros(n, n);
    tx_vinv_x_.zeros(m, m);
    pdv_.zeros(n, n, n_rands);
    logger_->info("Using {:d} samples", n);
}

double LinearMixedModel::ComputeLogLikelihood() const
{
    return -0.5
           * (logdet_v_ + log_det_sympd(tx_vinv_x_)
              + as_scalar(y_.t() * proj_y_));
}

dmat LinearMixedModel::ComputeZKZ(const sp_dmat& z, const dmat& k)
{
    bool z_identity = CheckIdentity(z);
    bool k_identity = CheckIdentity(k);

    if (k_identity)
    {
        return z_identity ? k : dmat(z * z.t());
    }
    return z_identity ? k : dmat(z * k * z.t());
}

dcube LinearMixedModel::ComputeZKZtR()
{
    auto n = Z_[0].n_rows;
    auto n_rands = covar_matrices_rand_.n_slices;
    dcube result(n, n, n_rands + 1, arma::fill::zeros);

    for (size_t i = 0; i < n_rands; ++i)
    {
        result.slice(i) = ComputeZKZ(Z_[i], covar_matrices_rand_.slice(i));
    }
    result.slice(n_rands).eye();

    return result;
}

void LinearMixedModel::ComputeV()
{
    v_.zeros();
    for (size_t i{0}; i < covar_matrices_rand_.n_slices; ++i)
    {
        v_ += sigma_.at(i) * pdv_.slice(i);
    }
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
    for (size_t i = 0; i < covar_matrices_rand_.n_slices; ++i)
    {
        pdv_.slice(i) = proj_ * zkztr_.slice(i);
    }
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
