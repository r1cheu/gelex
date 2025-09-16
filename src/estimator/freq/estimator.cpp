#include "gelex/estimator/freq/estimator.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <armadillo>

#include "gelex/model/freq/model.h"
#include "gelex/optim/optimizer.h"
#include "gelex/optim/policy.h"
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
    logger_.log_model_info(model, optimizer_name_, tol_, max_iter_);

    optimizer_.init(model);
    em_step(model, em_init);

    size_t iter{1};
    double time_cost{};
    logger_.log_iter_header(model);

    for (; iter <= max_iter_; ++iter)
    {
        {
            Timer timer{time_cost};
            if (optimizer_name_ == "AI")
            {
                optimizer_.step<AverageInformationPolicy>(model);
            }
        }

        logger_.log_iteration(
            iter, optimizer_.loglike_, model.effects(), time_cost);

        if (optimizer_.converged_)
        {
            converged_ = true;
            break;
        }
    }
    iter_count_ = iter;

    report_results(model, start);
}

void Estimator::em_step(GBLUP& model, bool em_init)
{
    if (em_init)
    {
        double time_cost{};
        {
            Timer timer{time_cost};
            optimizer_.step<ExpectationMaximizationPolicy>(model);
        }

        logger_.log_em_initialization(
            optimizer_.loglike_, model.effects(), time_cost);
    }
}

void Estimator::report_results(
    GBLUP& model,
    const std::chrono::steady_clock::time_point& start_time)
{
    compute_var_se(model, optimizer_);

    auto end = std::chrono::steady_clock::now();
    double elapsed_time
        = std::chrono::duration<double>(end - start_time).count();

    compute_beta(model);
    compute_u(model);

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
        = arma::diagvec(arma::sqrt(arma::inv_sympd(optimizer_.tx_vinv_x_)));
    logger_.log_fixed_effects(model, fixed_se);

    logger_.log_variance_components(model);

    auto [h2_se, sum_var] = compute_h2_se(model);
    logger_.log_heritability(model, h2_se, sum_var);
    logger_.log_results_footer();
}

void Estimator::compute_beta(GBLUP& model)
{
    model.fixed_->coeff = arma::inv_sympd(optimizer_.tx_vinv_x_)
                          * (model.fixed_->design_matrix.t() * optimizer_.v_
                             * model.phenotype_);
}

void Estimator::compute_u(GBLUP& model)
{
    for (auto& effect : model.random_)
    {
        effect.coeff
            = effect.design_matrix.t() * optimizer_.proj_y_ * effect.sigma;
    }

    for (auto& effect : model.genetic_)
    {
        effect.coeff = effect.genetic_relationship_matrix
                       * effect.design_matrix.t() * optimizer_.proj_y_
                       * effect.sigma;
    }
}

double Estimator::compute_aic(const GBLUP& model) const
{
    auto k
        = static_cast<double>(model.effects().size() + model.n_fixed_effects());
    double aic = (-2 * optimizer_.loglike_) + (2 * k);
    return aic;
}

double Estimator::compute_bic(const GBLUP& model) const
{
    auto k
        = static_cast<double>(model.effects().size() + model.n_fixed_effects());
    auto n = static_cast<double>(model.n_individuals_);
    double bic = (-2 * optimizer_.loglike_) + (k * std::log(n));
    return bic;
}

void Estimator::compute_var_se(GBLUP& model, const Optimizer& optim) const
{
    dvec var_se = arma::sqrt(arma::diagvec(-optim.hess_inv_));
    for (size_t i = 0; i < model.effects_.size(); ++i)
    {
        model.effects_[i].se = var_se.at(i);
    }
}

std::pair<std::vector<double>, double> Estimator::compute_h2_se(
    const GBLUP& model) const
{
    auto n = model.effects_.size();
    double sum_var = arma::sum(arma::dvec(model.effects_.values()));
    double sum_sq = sum_var * sum_var;

    std::vector<double> grad;
    std::vector<double> h2_se;
    h2_se.reserve(model.n_genetic_effects());

    for (const auto& g_effect : model.genetic_)
    {
        for (const auto& effect : model.effects_)
        {
            if (effect.name == g_effect.name)
            {
                grad.emplace_back((sum_var - effect.sigma) / sum_sq);
            }
            else
            {
                grad.emplace_back(-effect.sigma / sum_sq);
            }
        }
        dvec grad_vec(grad);
        grad.clear();
        h2_se.emplace_back(
            std::sqrt(
                arma::as_scalar(
                    grad_vec.t() * -optimizer_.hess_inv_ * grad_vec)));
    }
    return {h2_se, sum_var};
}

}  // namespace gelex
