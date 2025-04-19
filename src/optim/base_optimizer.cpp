#include "gelex/optim/base_optimizer.h"

#include <cstdint>
#include <cstdlib>

#include <fmt/color.h>
#include <fmt/ranges.h>
#include <armadillo>

#include "gelex/model/gblup.h"
#include "gelex/utils.h"

namespace gelex
{
void OptimizerBase::step(GBLUP& model)
{
    step_inner(model);
    check_convergence(model);
}

void OptimizerBase::init(GBLUP& model)
{
    phenotype_var_ = arma::var(model.phenotype());
    n_sigma_ = model.n_genetic_effects() + model.n_group_effects() + 1;
    n_individuals_ = model.n_individuals();

    if (model.sigma().empty())  // in case user refit the model or given a prior
                                // for sigma
    {
        sigma_ = dvec(
            n_sigma_,
            arma::fill::value(phenotype_var_ / static_cast<double>(n_sigma_)));
        model.set_sigma(sigma_);
    }
    else
    {
        sigma_ = model.sigma();
    }

    first_grad_ = arma::zeros<dvec>(n_sigma_);
    dvpy_ = arma::zeros<dmat>(n_individuals_, n_sigma_);
    proj_ = arma::zeros<dmat>(n_individuals_, n_individuals_);
    v_ = arma::zeros<dmat>(n_individuals_, n_individuals_);
    tx_vinv_x_
        = arma::zeros<dmat>(model.n_common_effects(), model.n_common_effects());
}

void OptimizerBase::prepare_proj(const GBLUP& model)
{
    compute_v(model);
    compute_proj(model);
}

void OptimizerBase::compute_v(const GBLUP& model)
{
    uint64_t sigma_idx{};
    v_.zeros();
    for (const sp_dmat& mat : model.group_cov_mats())
    {
        v_ += model.sigma().at(sigma_idx) * mat;
        ++sigma_idx;
    }
    for (const dmat& mat : model.genetic_cov_mats())
    {
        v_ += model.sigma().at(sigma_idx) * mat;
        ++sigma_idx;
    }
    v_.diag() += model.sigma().back();
}

void OptimizerBase::compute_proj(const GBLUP& model)
{
    const dmat& design_mat = model.design_mat_beta();
    logdet_v_ = v_inv_logdet(v_);
    dmat vinv_x = v_ * design_mat;
    tx_vinv_x_ = design_mat.t() * vinv_x;

    proj_
        = v_
          - vinv_x
                * solve(tx_vinv_x_, vinv_x.t(), arma::solve_opts::likely_sympd);
    proj_y_ = proj_ * model.phenotype();
}

double OptimizerBase::compute_loglike(const GBLUP& model)
{
    return -0.5
           * (logdet_v_ + log_det_sympd(tx_vinv_x_)
              + as_scalar(model.phenotype().t() * proj_y_));
}

void OptimizerBase::compute_dvpy(const GBLUP& model)
{
    uint64_t counts{};
    for (const sp_dmat& mat : model.group_cov_mats())
    {
        dvpy_.unsafe_col(counts) = (mat * proj_y_);
        ++counts;
    }

    for (const dmat& mat : model.genetic_cov_mats())
    {
        dvpy_.unsafe_col(counts) = (mat * proj_y_);
        ++counts;
    }
    dvpy_.unsafe_col(counts) = proj_y_;
}

void OptimizerBase::compute_first_grad(const GBLUP& model)
{
    compute_dvpy(model);
    uint64_t counts{};
    for (const sp_dmat& mat : model.group_cov_mats())
    {
        first_grad_.at(counts)
            = -0.5
              * (arma::trace(proj_ * mat)
                 - as_scalar(proj_y_.t() * dvpy_.unsafe_col(counts)));
        counts++;
    }

    for (const dmat& mat : model.genetic_cov_mats())
    {
        first_grad_.at(counts)
            = -0.5
              * (arma::trace(proj_ * mat)
                 - as_scalar(proj_y_.t() * dvpy_.unsafe_col(counts)));
        counts++;
    }
    first_grad_.back()
        = -0.5
          * (arma::trace(proj_)
             - arma::dot(dvpy_.unsafe_col(counts), dvpy_.unsafe_col(counts)));
}

double OptimizerBase::compute_sigma_diff(const GBLUP& model)
{
    double diff = norm(model.sigma() - sigma_) / norm(model.sigma());
    sigma_ = model.sigma();
    return diff;
}

double OptimizerBase::compute_loglike_diff(const GBLUP& model)
{
    double new_value = compute_loglike(model);
    double diff{new_value - loglike_};

    loglike_ = new_value;
    loglike_diff_ = diff;
    return diff;
}

void OptimizerBase::check_convergence(const GBLUP& model)
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

dvec OptimizerBase::constrain(const dvec& sigma, double y_var)
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
