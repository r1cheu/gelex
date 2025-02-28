#include "chenx/model/linear_mixed_model.h"

#include <stdexcept>
#include <vector>

#include <armadillo>

namespace chenx
{
LinearMixedModel::LinearMixedModel(
    dmat&& y,
    dmat&& X,
    dcube&& covar_matrices_rand,
    std::vector<std::string>&& random_effect_names)
    : y_{std::move(y)},
      X_{std::move(X)},
      zkzt_{covar_matrices_rand},
      random_effect_names_{std::move(random_effect_names)}
{
    y_var_ = arma::as_scalar(arma::cov(y_));  // currently only support scalar

    num_fixed_effects_ = X_.n_cols;
    num_individuals_ = X_.n_rows;
    num_random_effects_ = zkzt_.n_slices + 1;

    set_beta(dvec(num_fixed_effects_, arma::fill::zeros));
    pdv_.set_size(num_individuals_, num_individuals_, num_random_effects_);
    set_sigma(dvec(
        num_random_effects_,
        arma::fill::value(y_var_ / static_cast<double>(num_random_effects_))));
    random_effect_names_.emplace_back("e");
    U_.set_size(num_individuals_, num_random_effects_);
}

void LinearMixedModel::set_sigma(dvec&& sigma)
{
    sigma_ = std::move(sigma);
    ComputeV();
    ComputeProj();
    ComputePdV();
}

void LinearMixedModel::Reset()
{
    set_sigma(dvec(
        num_random_effects_,
        arma::fill::value(y_var_ / static_cast<double>(num_random_effects_))));
    set_beta(dvec(num_fixed_effects_, arma::fill::zeros));
};

double LinearMixedModel::ComputeLogLikelihood() const
{
    return -0.5
           * (logdet_v_ + log_det_sympd(tx_vinv_x_)
              + as_scalar(y_.t() * proj_y_));
}

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
