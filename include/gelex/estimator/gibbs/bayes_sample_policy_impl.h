#pragma once

#include "gelex/estimator/gibbs/bayes_sample_policy.h"

#include <armadillo>
#include <cmath>

#include "gelex/dist.h"
#include "gelex/estimator/gibbs/base.h"
#include "gelex/estimator/gibbs/bayes_b.h"
#include "gelex/estimator/gibbs/bayes_c.h"

#include "gelex/model/bayes_model.h"

namespace gelex
{

template <typename ModelType, typename Derived>
BaseGeneticSampler<ModelType, Derived>::BaseGeneticSampler(
    std::mt19937_64& random_generator,
    ModelType& model)
    : normal_{random_generator},
      genotype_mat_{model.genotype_mat()},
      cols_norm_{model.a_cols_norm()},
      cols_var_{model.a_cols_var()},
      gen_{random_generator}
{
}

GeneticSampler<BayesRR>::GeneticSampler(
    std::mt19937_64& random_generator,
    BayesRR& model)
    : BaseGeneticSampler<BayesRR, GeneticSampler<BayesRR>>(
          random_generator,
          model),
      chisq_{
          random_generator,
          model.priors().sigma_g().nu,
          static_cast<double>(model.a().n_elem),
          model.priors().sigma_g().s2}
{
}

void GeneticSampler<BayesRR>::SampleImpl(BayesRR& model, dvec& y_adj)
{
    dvec& a = model.a();  // NOLINT
    double sigma_a = model.sigma_a();
    const double sigma_e = model.sigma_e();
    const dvec& inv_scaler = 1.0 / (cols_norm_ + sigma_e / sigma_a);

    for (uint64_t i = 0; i < a.n_elem; ++i)
    {
        if (cols_var_.at(i) == 0.0)
        {
            continue;
        }
        const double old_i = a.at(i);
        const dvec& col_i = genotype_mat_.col(i);
        const double col_norm = cols_norm_.at(i);
        const double inv_scaler_i = inv_scaler.at(i);

        double rhs = ComputeRhs(col_i, y_adj, old_i, col_norm);
        double new_i
            = normal_(rhs * inv_scaler_i, sqrt(sigma_e * inv_scaler_i));

        a.at(i) = new_i;
        daxpy_auto(y_adj, col_i, old_i - new_i);
    }
    model.set_sigma_a(chisq_(arma::dot(a, a)));
}

GeneticSampler<BayesA>::GeneticSampler(
    std::mt19937_64& random_generator,
    BayesA& model)
    : BaseGeneticSampler<BayesA, GeneticSampler<BayesA>>(
          random_generator,
          model),
      chisq_{
          random_generator,
          model.priors().sigma_g().nu,
          1,
          model.priors().sigma_g().s2}

{
}

void GeneticSampler<BayesA>::SampleImpl(BayesA& model, dvec& y_adj)
{
    dvec& a = model.a();  // NOLINT
    dvec& sigma_a = model.sigma_a();
    const double sigma_e = model.sigma_e();

    for (uint64_t i = 0; i < a.n_elem; ++i)
    {
        if (cols_var_.at(i) == 0.0)
        {
            continue;
        }
        const double old_i = a.at(i);
        const dvec& col_i = genotype_mat_.col(i);
        const double col_norm = cols_norm_.at(i);
        double inv_scaler = 1.0 / (col_norm + sigma_e / sigma_a.at(i));

        double rhs = ComputeRhs(col_i, y_adj, old_i, col_norm);
        double new_i = normal_(rhs * inv_scaler, sqrt(sigma_e * inv_scaler));

        a.at(i) = new_i;
        daxpy_auto(y_adj, col_i, old_i - new_i);
        sigma_a.at(i) = chisq_(new_i * new_i);
    }
}

GeneticSampler<BayesB>::GeneticSampler(
    std::mt19937_64& random_generator,
    BayesB& model)
    : BaseGeneticSampler<BayesB, GeneticSampler<BayesB>>(
          random_generator,
          model),
      uniform_{random_generator},
      chisq_{
          random_generator,
          model.priors().sigma_g().nu,
          1,
          model.priors().sigma_g().s2},
      snp_tracker{model.a().n_elem, arma::fill::zeros},
      fold_(model.pi().n_elem, arma::fill::zeros)
{
}

void GeneticSampler<BayesB>::SampleImpl(BayesB& model, dvec& y_adj)
{
    BayesBKernel<BayesB, false>(
        model,
        y_adj,
        genotype_mat_,
        cols_norm_,
        cols_var_,
        normal_,
        chisq_,
        uniform_,
        snp_tracker,
        fold_,
        gen_);
}

GeneticSampler<BayesBpi>::GeneticSampler(
    std::mt19937_64& random_generator,
    BayesBpi& model)
    : BaseGeneticSampler<BayesBpi, GeneticSampler<BayesBpi>>(
          random_generator,
          model),
      uniform_{random_generator},
      chisq_{
          random_generator,
          model.priors().sigma_g().nu,
          1,
          model.priors().sigma_g().s2},
      snp_tracker{model.a().n_elem, arma::fill::zeros},
      fold_(model.pi().n_elem, arma::fill::zeros)
{
}

void GeneticSampler<BayesBpi>::SampleImpl(BayesBpi& model, dvec& y_adj)
{
    BayesBKernel<BayesBpi, true>(
        model,
        y_adj,
        genotype_mat_,
        cols_norm_,
        cols_var_,
        normal_,
        chisq_,
        uniform_,
        snp_tracker,
        fold_,
        gen_);
}

GeneticSampler<BayesC>::GeneticSampler(
    std::mt19937_64& random_generator,
    BayesC& model)
    : BaseGeneticSampler<BayesC, GeneticSampler<BayesC>>(
          random_generator,
          model),
      uniform_{random_generator},
      chisq_{
          random_generator,
          model.priors().sigma_g().nu,
          1,
          model.priors().sigma_g().s2},
      snp_tracker{model.a().n_elem, arma::fill::zeros},
      fold_(model.pi().n_elem, arma::fill::zeros)
{
}

void GeneticSampler<BayesC>::SampleImpl(BayesC& model, dvec& y_adj)
{
    BayesCKernel<BayesC, false>(
        model,
        y_adj,
        genotype_mat_,
        cols_norm_,
        cols_var_,
        normal_,
        chisq_,
        uniform_,
        snp_tracker,
        fold_,
        gen_);
}

GeneticSampler<BayesCpi>::GeneticSampler(
    std::mt19937_64& random_generator,
    BayesCpi& model)
    : BaseGeneticSampler<BayesCpi, GeneticSampler<BayesCpi>>(
          random_generator,
          model),
      uniform_{random_generator},
      chisq_{
          random_generator,
          model.priors().sigma_g().nu,
          1,
          model.priors().sigma_g().s2},
      snp_tracker{model.a().n_elem, arma::fill::zeros},
      fold_(model.pi().n_elem, arma::fill::zeros)
{
}

void GeneticSampler<BayesCpi>::SampleImpl(BayesCpi& model, dvec& y_adj)
{
    BayesCKernel<BayesCpi, true>(
        model,
        y_adj,
        genotype_mat_,
        cols_norm_,
        cols_var_,
        normal_,
        chisq_,
        uniform_,
        snp_tracker,
        fold_,
        gen_);
}
}  // namespace gelex
