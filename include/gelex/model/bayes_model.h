#pragma once
#include <armadillo>
#include <optional>
#include <string>

#include "gelex/model/bayes_model_policy.h"
#include "gelex/model/bayes_prior.h"

namespace gelex
{
using arma::dmat;
using arma::dvec;
using arma::uvec;

template <GeneticPolicy Genetic>
class BayesianModel
{
   public:
    BayesianModel(
        dvec phenotype,
        dmat genotype_mat,
        std::optional<dmat> design_mat_beta,
        std::optional<dmat> design_mat_r)
        : phenotype_{std::move(phenotype)},
          genotype_mat_{std::move(genotype_mat)},
          design_mat_beta_{std::move(design_mat_beta)},
          design_mat_r_{std::move(design_mat_r)},
          priors_{Genetic::init_pi()}
    {
        a_.zeros(genotype_mat_.n_cols);
        a_cols_norm_ = sum_square(genotype_mat_);
        a_cols_var_ = cols_var(genotype_mat_);
        n_var_0_ = arma::sum(a_cols_var_ == 0.0);

        sigma_a_ = Genetic::init_sigma(genotype_mat_.n_cols);  // NOLINT

        if (design_mat_beta_)
        {
            beta_.zeros(design_mat_beta_->n_cols);
            beta_cols_norm_ = sum_square(*design_mat_beta_);
        }

        if (design_mat_r_)
        {
            r_.zeros(design_mat_r_->n_cols);
            r_cols_norm_ = sum_square(*design_mat_r_);
        }
    };

    static constexpr std::string name = Genetic::name;
    static constexpr bool has_pi = Genetic::has_pi;
    static constexpr bool fixed_pi = Genetic::fixed_pi;
    using sigma_t = Genetic::sigma_t;
    const dvec& phenotype() const { return phenotype_; }
    const dmat& genotype_mat() const { return genotype_mat_; }
    const std::optional<dmat>& design_mat_beta() const
    {
        return design_mat_beta_;
    }
    const std::optional<dmat>& design_mat_r() const { return design_mat_r_; }
    Priors& priors() { return priors_; }
    const Priors& priors() const { return priors_; }

    const dvec& pi() const { return pi_; }
    void set_pi(dvec new_pi) { pi_ = std::move(new_pi); }

    double mu() const { return mu_; }
    void set_mu(double mu) { mu_ = mu; }  // NOLINT

    dvec& a() { return a_; }
    dvec& beta() { return beta_; }
    dvec& r() { return r_; }

    template <typename T = sigma_t>
    auto sigma_a() const ->
        typename std::enable_if_t<std::is_same_v<T, double>, double>
    {
        return sigma_a_;
    }

    template <typename T = sigma_t>
    auto sigma_a() -> typename std::enable_if_t<!std::is_same_v<T, double>, T&>
    {
        return sigma_a_;
    }

    template <typename T = sigma_t>
    auto set_sigma_a(double sigma_a) const ->
        typename std::enable_if_t<std::is_same_v<T, double>, void>
    {
        sigma_a_ = sigma_a;
    }

    double sigma_r() { return sigma_r_; }
    void set_sigma_r(double sigma_r) { sigma_r_ = sigma_r; }
    double sigma_e() { return sigme_e_; }
    void set_sigma_e(double sigma_e) { sigme_e_ = sigma_e; }

    const dvec& a_cols_var() const { return a_cols_var_; }
    const dvec& a_cols_norm() const { return a_cols_norm_; }
    uint64_t n_var_0() const { return n_var_0_; }

    const dvec& beta_cols_norm() const { return beta_cols_norm_; }
    const dvec& r_cols_norm() const { return r_cols_norm_; }

    bool has_env() const { return design_mat_r_.has_value(); }
    bool has_beta() const { return design_mat_beta_.has_value(); }

    void set_priors() { pi_ = priors_.pi(); }

   protected:
    static dvec sum_square(const dmat& mat)
    {
        dvec out(mat.n_cols);
#pragma omp parallel for default(none) shared(mat, out)
        for (size_t i = 0; i < mat.n_cols; ++i)
        {
            out.at(i) = arma::dot(mat.unsafe_col(i), mat.unsafe_col(i));
        }
        return out;
    }
    static dvec cols_var(const dmat& mat)
    {
        dvec out(mat.n_cols);
#pragma omp parallel for default(none) shared(mat, out)
        for (size_t i = 0; i < mat.n_cols; ++i)
        {
            out.at(i) = arma::var(mat.unsafe_col(i));
        }
        return out;
    }

   private:
    dvec phenotype_;     // Phenotype vector
    dmat genotype_mat_;  // Genotype matrix
    std::optional<dmat>
        design_mat_beta_;               // Design matrix for beta coefficients
    std::optional<dmat> design_mat_r_;  // Design matrix for rironmental factors
    Priors priors_;

    dvec pi_;  // proportion of each group

    double mu_{};
    dvec a_;
    dvec beta_;
    dvec r_;

    dvec a_cols_var_;     // variance of each column in genotype matrix
    uint64_t n_var_0_{};  // number of columns with zero var
    dvec a_cols_norm_;

    dvec beta_cols_norm_;
    dvec r_cols_norm_;

    sigma_t sigma_a_;
    double sigma_r_{};
    double sigme_e_{};
};

using BayesA = BayesianModel<BayesAPolicy>;
using BayesB = BayesianModel<BayesBPolicy>;
using BayesBpi = BayesianModel<BayesBpiPolicy>;
using BayesC = BayesianModel<BayesCPolicy>;
using BayesCpi = BayesianModel<BayesCpiPolicy>;
using BayesR = BayesianModel<BayesRPolicy>;
using BayesRR = BayesianModel<BayesRRPolicy>;
}  // namespace gelex
