#include "gelex/model/bayes.h"
#include <cmath>
#include "armadillo"

namespace gelex
{

BayesRR::BayesRR(
    dvec&& phenotype,
    dmat&& design_matrix_beta,
    dmat&& genotype_matrix,
    uint64_t seed,
    uint64_t n_iter,
    uint64_t n_burn)
    : phenotype_(std::move(phenotype)),
      design_matrix_beta_(std::move(design_matrix_beta)),
      genotype_matrix_(std::move(genotype_matrix)),
      n_iter_{n_iter},
      n_burn_{n_burn}

{
    rng_.seed(seed);
}

void BayesRR::set_hyperparameters(double nu_e, double nu_a)
{
    nu_e_ = nu_e;
    nu_a_ = nu_a;
    int n = 0;
}

void BayesRR::SamplingBeta(uint64_t index)
{
    double col_norm2 = b_cols_norm2_.at(index);
    double new_beta = SampleElement(
        design_matrix_beta_.col(index),
        col_norm2,
        beta_.at(index),
        sigma_e_,
        0.0);  // No regularization ratio for beta
    rss_ += arma::as_scalar(
        genotype_matrix_.col(index).t() * (new_beta - beta_.at(index)));
    beta_.at(index) = new_beta;
}

void BayesRR::SamplingSnpEffect(uint64_t index)
{
    double col_norm2 = g_cols_norm2_.at(index);
    double sigme_ratio = sigma_e_ / sigma_a_;
    double new_snp_eff = SampleElement(
        genotype_matrix_.col(index),
        col_norm2,
        snp_effect_.at(index),
        sigma_e_,
        sigme_ratio);
    rss_ += arma::as_scalar(
        genotype_matrix_.col(index).t()
        * (new_snp_eff - snp_effect_.at(index)));
    snp_effect_.at(index) = new_snp_eff;
}

double BayesRR::SampleElement(
    const arma::vec& col_vector,
    double col_norm2,
    double current_effect,
    double sigma,
    double regularization_ratio)
{
    double denominator = col_norm2 + regularization_ratio;
    double mean = arma::as_scalar(
        (col_vector.t() * rss_ + col_norm2 * current_effect) / denominator);
    double std_dev = std::sqrt(sigma / denominator);
    return mean + (std_dev * normal_(rng_));
}

void BayesRR::init_params()
{
    mu_ = arma::mean(phenotype_);
    sigma_e_ = arma::var(phenotype_);
    sigma_a_ = sigma_e_ / 10.0;

    beta_ = arma::zeros(design_matrix_beta_.n_cols);
    snp_effect_ = arma::zeros(genotype_matrix_.n_cols);

    g_cols_norm2_ = arma::vecnorm(genotype_matrix_, 2, 0);
    b_cols_norm2_ = arma::vecnorm(design_matrix_beta_, 2, 0);

    s_a0_ = 1.0;
}
}  // namespace gelex
