#pragma once

#include <armadillo>

#include "gelex/model/gblup.h"
#include "gelex/utils.h"

namespace gelex
{
class OptimizerBase
{
    static constexpr double DEFAULT_TOL = 1e-8;

   public:
    explicit OptimizerBase(double tol = DEFAULT_TOL)
        : tol_{tol}, logger_{Logger::logger()}
    {
    }
    OptimizerBase(const OptimizerBase&) = default;
    OptimizerBase(OptimizerBase&&) noexcept = default;
    OptimizerBase& operator=(const OptimizerBase&) = default;
    OptimizerBase& operator=(OptimizerBase&&) noexcept = default;
    virtual ~OptimizerBase() = default;

    double tol() const noexcept { return tol_; }
    bool converged() const noexcept { return converged_; }
    double loglike() const noexcept { return loglike_; }

    void init(GBLUP& model);
    void set_tol(double tol) noexcept { tol_ = tol; }
    void step(GBLUP& model);

    double phenotype_var() const noexcept { return phenotype_var_; }

    const dmat& v() const noexcept { return v_; }
    const dmat& tx_vinv_x() const noexcept { return tx_vinv_x_; }
    const dmat& proj_y() const noexcept { return proj_y_; }

   protected:
    dvec constrain(const dvec& sigma, double y_var);
    void prepare_proj(const GBLUP& model);

    const dmat& dvpy() const noexcept { return dvpy_; }
    const dmat& proj() const noexcept { return proj_; }
    const dvec& first_grad() const noexcept { return first_grad_; }
    double loglike_diff() const noexcept { return loglike_diff_; }

    void compute_dvpy(const GBLUP& model);
    void compute_pdv(const GBLUP& model);
    void compute_first_grad(const GBLUP& model);

    void check_convergence(const GBLUP& model);

   private:
    void compute_v(const GBLUP& model);
    void compute_proj(const GBLUP& model);
    virtual void step_inner(GBLUP& model) = 0;

    double compute_loglike(const GBLUP& model);
    double compute_sigma_diff(const GBLUP& model);
    double compute_loglike_diff(const GBLUP& model);

    double tol_;
    double phenotype_var_{};

    bool converged_{};
    dvec old_sigma_;
    double loglike_{};
    double loglike_diff_{};
    double logdet_v_{};

    dmat dvpy_;

    dvec first_grad_;
    dvec proj_y_;
    dmat v_, proj_, tx_vinv_x_;

    std::shared_ptr<spdlog::logger> logger_;
};

double v_inv_logdet(dmat&);
void solve_sympd(dmat& A, dmat& B);
}  // namespace gelex
