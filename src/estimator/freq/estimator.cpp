#include "gelex/estimator/freq/estimator.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <armadillo>

#include "gelex/model/freq/effects.h"
#include "gelex/model/freq/model.h"
#include "gelex/optim/base.h"
#include "gelex/optim/optimizers.h"
#include "gelex/utils/formatter.h"
#include "gelex/utils/utils.h"

namespace gelex
{

using arma::dmat;
using arma::dvec;

Estimator::Estimator(std::string_view optimizer, size_t max_iter, double tol)
    : max_iter_{max_iter}
{
    set_optimizer(optimizer, tol);
}

void Estimator::set_optimizer(std::string_view optimizer, double tol)
{
    std::string opt_lower = ToLowercase(optimizer);
    if (opt_lower == "ai")
    {
        optimizer_name_ = "AI";
        tol_ = tol;
    }
    else
    {
        throw std::invalid_argument(
            "Unknown optimizer: " + std::string(optimizer)
            + ", AI(Average Information) is supported.");
    }
}

void Estimator::fit(GBLUP& model, bool em_init, bool verbose)
{
    auto start = std::chrono::steady_clock::now();

    logger_.set_verbose(verbose);
    logger_.log_model_information(model, optimizer_name_, tol_, max_iter_);

    initialize_optimizer(model, em_init);
    run_optimization_loop(model);
    report_results(model, start);
}

void Estimator::initialize_optimizer(GBLUP& model, bool em_init)
{
    double time_cost{};
    if (em_init)
    {
        ExpectationMaximizationOptimizer em_optimizer{tol_};
        em_optimizer.init(model);

        {
            Timer timer{time_cost};
            em_optimizer.step(model);
        }

        logger_.log_em_initialization(
            em_optimizer.loglike(), model.random(), time_cost);

        optimizer_ = std::make_unique<AverageInformationOptimizer>(
            std::move(static_cast<OptimizerBase&>(em_optimizer)));
    }
    else
    {
        optimizer_ = std::make_unique<AverageInformationOptimizer>(tol_);
        optimizer_->init(model);
    }
}

void Estimator::run_optimization_loop(GBLUP& model)
{
    size_t iter{1};
    double time_cost{};

    for (; iter <= max_iter_; ++iter)
    {
        {
            Timer timer{time_cost};
            optimizer_->step(model);
        }

        logger_.log_iteration(
            iter, optimizer_->loglike(), model.random(), time_cost);

        if (optimizer_->converged())
        {
            converged_ = true;
            break;
        }
    }
    iter_count_ = iter;
}

void Estimator::report_results(
    GBLUP& model,
    const std::chrono::steady_clock::time_point& start_time)
{
    model.random().set_se(compute_se(model.random().hess_inv()));
    auto end = std::chrono::steady_clock::now();
    double elapsed_time
        = std::chrono::duration<double>(end - start_time).count();

    model.fixed().beta = compute_beta(model);

    logger_.log_results_header();
    logger_.log_convergence_status(
        converged_,
        iter_count_,
        max_iter_,
        elapsed_time,
        compute_aic(model),
        compute_bic(model));

    // Compute fixed effects standard errors
    dvec fixed_se
        = arma::diagvec(arma::sqrt(arma::inv_sympd(optimizer_->tx_vinv_x())));
    logger_.log_fixed_effects(model, fixed_se);

    logger_.log_variance_components(model);

    auto [h2_se, sum_var] = compute_h2_se(model.random());
    logger_.log_heritability(model, h2_se, sum_var);

    logger_.log_results_footer();
    compute_u(model);
}

dvec Estimator::compute_beta(GBLUP& model)
{
    auto compute_beta_visitor = [this, &model](const auto& mat) -> dvec
    {
        return arma::inv_sympd(optimizer_->tx_vinv_x())
               * (mat.t() * optimizer_->v() * model.phenotype());
    };
    return std::visit(compute_beta_visitor, model.fixed().design_mat);
}

dmat Estimator::compute_u(GBLUP& model)
{
    dmat U{
        model.n_individuals(),
        model.random().sigma().n_elem,
        arma::fill::zeros};

    size_t idx{};
    for (const auto& effect : model.random())
    {
        std::visit(
            [&](const auto& cov)
            { U.unsafe_col(idx) = cov * optimizer_->proj_y() * effect.sigma; },
            effect.cov_mat);
        ++idx;
    }
    return U;
}

double Estimator::compute_aic(GBLUP& model)
{
    auto k
        = static_cast<double>(model.random().size() + model.n_fixed_effects());
    double aic = (-2 * optimizer_->loglike()) + (2 * k);
    return aic;
}

double Estimator::compute_bic(GBLUP& model)
{
    auto k
        = static_cast<double>(model.random().size() + model.n_fixed_effects());
    auto n = static_cast<double>(model.n_individuals());
    double bic = (-2 * optimizer_->loglike()) + (k * std::log(n));
    return bic;
}

dvec compute_se(const dmat& hess_inv)
{
    return arma::sqrt(arma::diagvec(-hess_inv));
}

std::pair<std::vector<double>, double> compute_h2_se(
    const RandomEffectManager& effects)
{
    auto n = effects.size();
    double sum_var = arma::sum(effects.sigma());
    double sum_sq = sum_var * sum_var;

    dvec grad = arma::zeros<dvec>(n);
    std::vector<double> h2_se;
    h2_se.reserve(effects.n_genetic_effects());
    for (auto i : effects.genetic_indices())
    {
        for (size_t j{}; j < n; ++j)
        {
            if (i == j)
            {
                grad[j] = (sum_var - effects[i].sigma) / sum_sq;
            }
            else
            {
                grad[j] = -effects[i].sigma / sum_sq;
            }
        }
        h2_se.emplace_back(
            std::sqrt(arma::as_scalar(grad.t() * -effects.hess_inv() * grad)));
    }
    return {h2_se, sum_var};
}

}  // namespace gelex
