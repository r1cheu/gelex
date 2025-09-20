#pragma once

#include <armadillo>

#include "gelex/logger.h"
#include "gelex/model/freq/model.h"

namespace gelex
{
class Optimizer
{
    friend class ExpectationMaximizationPolicy;
    friend class AverageInformationPolicy;
    friend class Estimator;

   public:
    explicit Optimizer(double tol = 1e-8)
        : tol_{tol}, logger_{gelex::logging::get()}
    {
    }
    void init(GBLUP& model);
    void set_tol(double tol) { tol_ = tol; }

    template <typename OptimPolicy>
    void step(GBLUP& model)
    {
        arma::dvec new_sigma = OptimPolicy::apply(*this, model);
        model.set_var_comp(new_sigma);
        check_convergence(model);
    }

   private:
    void compute_proj(const GBLUP& model);
    void compute_dvpy(const GBLUP& model);
    void compute_pdv(const GBLUP& model);
    void compute_first_grad(const GBLUP& model);

    arma::dvec constrain(const arma::dvec& sigma, double y_var);
    void check_convergence(const GBLUP& model);

    double compute_loglike(const GBLUP& model);
    double compute_sigma_diff(const GBLUP& model);
    double compute_loglike_diff(const GBLUP& model);

    template <typename EffectVisitor>
    void visitor_effects(
        const GBLUP& model,
        EffectVisitor visitor,
        size_t start = 1) const
    {
        size_t idx = start;
        for (const auto& eff : model.random_)
        {
            visitor(eff, idx++);
        }
        for (const auto& eff : model.genetic_)
        {
            visitor(eff, idx++);
        }
        for (const auto& eff : model.gxe_)
        {
            visitor(eff, idx++);
        }
    }

    double tol_;
    double phenotype_var_{};

    bool converged_{};
    arma::dvec old_sigma_;
    double loglike_{};
    double loglike_diff_{};
    double logdet_v_{};

    arma::dmat dvpy_;

    arma::dvec first_grad_;
    arma::dvec proj_y_;
    arma::dmat v_, proj_, tx_vinv_x_;

    arma::dmat hess_inv_;

    std::shared_ptr<spdlog::logger> logger_;
};

double v_inv_logdet(arma::dmat&);
void solve_sympd(arma::dmat& A, arma::dmat& B);
}  // namespace gelex
