#include "chenx/model/linear_mixed_model.h"

#include "chenx/utils.h"

namespace chenx
{
using namespace arma;
LinearMixedModel::LinearMixedModel(
    dvec&& y,
    dmat&& X,
    const uvec& z_index,
    dcube&& rands,
    std::vector<std::string>&& rand_names)
    : y_{std::move(y)},
      X_{std::move(X)},
      rands_{std::move(rands)},
      rand_names_{std::move(rand_names)}
{
    uword n{X.n_rows}, k{X.n_cols},
        n_rands{rands.n_slices + 1};  // add 1 for sigma_e
    y_var_ = var(y);
    Z_ = CreateZ(n_rands - 1, z_index, rands_.n_rows);
    zkztr_ = ComputeZKZtR(Z_, rands);
    beta_.zeros(k);
    sigma_.set_size(n_rands).fill(var(y) / static_cast<double>(n_rands));

    proj_y_.zeros(n);
    v_.zeros(n, n);
    txvx_.zeros(k, k);
    pdv_.zeros(n, n, n_rands);
}

double LinearMixedModel::ComputeLogLikelihood() const
{
    return -0.5
           * (logdet_v_ + log_det_sympd(txvx_) + as_scalar(y_.t() * proj_y_));
}

std::vector<sp_dmat> LinearMixedModel::CreateZ(
    const uword& n_z,
    const Col<uword>& z_index,
    const uword& n)
{
    std::vector<sp_dmat> z;
    z.reserve(n_z);
    for (size_t i = 0; i < n_z; ++i)
    {
        z.emplace_back(speye<sp_dmat>(n, n).cols(z_index).t());
    }
    return z;
}

std::vector<sp_dmat> LinearMixedModel::CreateZ(const uword& n_z, const uword& n)
{
    std::vector<sp_dmat> z;
    z.reserve(n_z);
    for (size_t i = 0; i < n_z; ++i)
    {
        z.emplace_back(speye<sp_dmat>(n, n));
    }
    return z;
}

dmat LinearMixedModel::ComputeZKZ(const sp_dmat& z, const dmat& k)
{
    bool z_identity = check_identity(z);
    bool k_identity = check_identity(k);

    if (k_identity)
    {
        return z_identity ? k : dmat(z * z.t());
    }
    else
    {
        return z_identity ? k : dmat(z * k * z.t());
    }
}

dcube LinearMixedModel::ComputeZKZtR(
    const std::vector<sp_dmat>& z,
    const dcube& k)
{
    auto n = z[0].n_rows;
    dcube result(n, n, k.n_slices + 1, fill::zeros);

    for (size_t i = 0; i < k.n_slices; ++i)
    {
        result.slice(i) = ComputeZKZ(z[i], k.slice(i));
    }

    result.slice(k.n_slices).eye();
    return result;
}

void LinearMixedModel::ComputeV()
{
    v_.zeros();
    for (size_t i{0}; i < rands_.n_slices; ++i)
    {
        v_ += sigma_.at(i) * pdv_.slice(i);
    }
}

void LinearMixedModel::ComputeProj()
{
    logdet_v_ = VinvLogdet(v_);
    dmat vx = v_ * X_;
    txvx_ = X_.t() * vx;
    proj_ = v_ - vx * solve(txvx_, vx.t(), solve_opts::likely_sympd);
    proj_y_ = proj_ * y_;
}

void LinearMixedModel::ComputePdV()
{
    for (size_t i = 0; i < rands_.n_slices; ++i)
    {
        pdv_.slice(i) = proj_ * zkztr_.slice(i);
    }
}

double LinearMixedModel::VinvLogdet(dmat& V)
{
    char uplo = 'L';
    int n = V.n_cols;
    int info;
    lapack::potrf(&uplo, &n, V.memptr(), &n, &info);
    if (info != 0)
    {
        throw std::runtime_error("V Matrix is not symmetric positive definite");
    }
    double logdet = accu(2.0 * log(diagvec(V)));
    lapack::potri(&uplo, &n, V.memptr(), &n, &info);

    if (info != 0)
    {
        throw std::runtime_error("Error during inverse V matrix");
    }
    V = arma::symmatl(V);
    return logdet;
}

}  // namespace chenx
