#include "gelex/optim/optimizer.h"

#include <fmt/color.h>
#include <fmt/ranges.h>
#include <armadillo>

#include "gelex/model/freq/model.h"

namespace gelex
{
using arma::dmat;
using arma::dvec;

void Optimizer::init(GBLUP& model)
{
    phenotype_var_ = arma::var(model.phenotype_);
    model.init_var_comp(phenotype_var_);
    old_sigma_ = model.var_comp();

    size_t n_sigma = old_sigma_.n_elem;
    auto n_individuals = model.n_individuals_;

    first_grad_ = arma::zeros<dvec>(n_sigma);
    dvpy_ = arma::zeros<dmat>(n_individuals, n_sigma);
    proj_ = arma::zeros<dmat>(n_individuals, n_individuals);
    v_ = arma::zeros<dmat>(n_individuals, n_individuals);
    tx_vinv_x_
        = arma::zeros<dmat>(model.n_fixed_effects(), model.n_fixed_effects());
}

void Optimizer::compute_proj(const GBLUP& model)
{
    v_.zeros();
    v_.diag() += model.residual_.sigma;
    auto accumulate = [this](const auto& effect, size_t place_holder)
    { v_ += effect.covariance_matrix * effect.sigma; };
    visitor_effects(model, accumulate);

    logdet_v_ = v_inv_logdet(v_);
    dmat vinv_x = v_ * model.fixed_->design_matrix;
    tx_vinv_x_ = model.fixed_->design_matrix.t() * vinv_x;

    proj_
        = v_
          - vinv_x
                * solve(tx_vinv_x_, vinv_x.t(), arma::solve_opts::likely_sympd);
    proj_y_ = proj_ * model.phenotype();
}

double Optimizer::compute_loglike(const GBLUP& model)
{
    return -0.5
           * (logdet_v_ + log_det_sympd(tx_vinv_x_)
              + as_scalar(model.phenotype_.t() * proj_y_));
}

void Optimizer::compute_dvpy(const GBLUP& model)
{
    dvpy_.unsafe_col(0) = proj_y_;
    auto apply_effects = [this](const auto& effect, size_t idx)
    { dvpy_.unsafe_col(idx) = effect.covariance_matrix * proj_y_; };
    visitor_effects(model, apply_effects);
}

void Optimizer::compute_first_grad(const GBLUP& model)
{
    compute_dvpy(model);
    first_grad_.at(0)
        = -0.5
          * (arma::trace(proj_) - as_scalar(proj_y_.t() * dvpy_.unsafe_col(0)));
    auto compute_first_grad = [this](const auto& effect, size_t idx)
    {
        first_grad_.at(idx)
            = -0.5
              * (arma::trace(proj_ * effect.covariance_matrix)
                 - as_scalar(proj_y_.t() * dvpy_.unsafe_col(idx)));
    };
    visitor_effects(model, compute_first_grad);
}

double Optimizer::compute_sigma_diff(const GBLUP& model)
{
    dvec new_sigma = model.var_comp();
    double diff = norm(new_sigma - old_sigma_) / norm(new_sigma);
    old_sigma_ = new_sigma;
    return diff;
}

double Optimizer::compute_loglike_diff(const GBLUP& model)
{
    double new_value = compute_loglike(model);
    double diff{new_value - loglike_};

    loglike_ = new_value;
    loglike_diff_ = diff;
    return diff;
}

void Optimizer::check_convergence(const GBLUP& model)
{
    double sigma_diff = compute_sigma_diff(model);
    double loglike_diff = compute_loglike_diff(model);
    bool negative = loglike_diff < 0.0;
    loglike_diff = std::abs(loglike_diff);

    if (sigma_diff < tol_
        && (loglike_diff < 1e-4 || (negative && loglike_diff < 1e-2)))
    {
        converged_ = true;
    }
}

dvec Optimizer::constrain(const dvec& sigma, double y_var)
{
    constexpr double constr_scale = 1e-6;
    arma::vec constrained_sigma = sigma;
    arma::uvec constrained = arma::find(sigma < 0.0);
    if (constrained.is_empty())
    {
        return constrained_sigma;
    }
    arma::uvec unconstrained = arma::find(sigma > 0.0);
    constrained_sigma(constrained).fill(y_var * constr_scale);

    double delta = arma::accu(y_var * constr_scale - sigma(constrained))
                   / static_cast<double>(unconstrained.n_elem);

    for (const auto& i : unconstrained)
    {
        if (sigma(i) > delta)
        {
            constrained_sigma.at(i) -= delta;
        }
    }

    if (constrained.n_elem > sigma.n_elem / 2)
    {
        logger_->warn(
            "Half of the variance components are constrained! The "
            "estimate is not reliable.");
    }

    return constrained_sigma;
}

double v_inv_logdet(dmat& V)
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

void solve_sympd(dmat& A, dmat& B)  // NOLINT
{
    char uplo = 'L';
    int n = static_cast<int>(A.n_cols);
    int k = static_cast<int>(B.n_cols);
    int info = 0;
    int lda = n;
    int ldb = n;
    arma::dposv_(&uplo, &n, &k, A.memptr(), &lda, B.memptr(), &ldb, &info);
}
}  // namespace gelex
