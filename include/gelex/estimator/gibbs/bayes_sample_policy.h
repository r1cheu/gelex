#pragma once

#include <armadillo>
#include <cmath>

#include "gelex/dist.h"
#include "gelex/estimator/gibbs/base.h"
#include "gelex/model/bayes_model.h"

namespace gelex
{

using dvec = arma::vec;

template <typename ModelType, typename Derived>
class BaseGeneticSampler
{
   private:
    Normal normal_;
    const dmat& genotype_mat_;
    const dvec& cols_norm_;
    const dvec& cols_var_;
    std::mt19937_64& gen_;
    BaseGeneticSampler(std::mt19937_64& random_generator, ModelType& model);

   public:
    BaseGeneticSampler(const BaseGeneticSampler&) = delete;  // NOLINT
    BaseGeneticSampler(BaseGeneticSampler&&) = delete;       // NOLINT
    BaseGeneticSampler& operator=(const BaseGeneticSampler&) = delete;
    BaseGeneticSampler& operator=(BaseGeneticSampler&&) = delete;
    ~BaseGeneticSampler() = default;

    void sample(ModelType& model, dvec& y_adj)
    {
        static_cast<Derived*>(this)->sample_impl(model, y_adj);
    }

    friend Derived;
};

template <typename ModelType>
class GeneticSampler
{
};

template <>
class GeneticSampler<BayesRR>
    : public BaseGeneticSampler<BayesRR, GeneticSampler<BayesRR>>
{
   public:
    GeneticSampler(std::mt19937_64& random_generator, BayesRR& model);
    void sample_impl(BayesRR& model, dvec& y_adj);

   private:
    ScaleInvChiSq chisq_;
};

template <>
class GeneticSampler<BayesA>
    : public BaseGeneticSampler<BayesA, GeneticSampler<BayesA>>
{
   public:
    GeneticSampler(std::mt19937_64& random_generator, BayesA& model);
    void sample_impl(BayesA& model, dvec& y_adj);

   private:
    ScaleInvChiSq chisq_;
};

template <>
class GeneticSampler<BayesB>
    : public BaseGeneticSampler<BayesB, GeneticSampler<BayesB>>
{
   public:
    GeneticSampler(std::mt19937_64& random_generator, BayesB& model);

    void sample_impl(BayesB& model, dvec& y_adj);

   private:
    ScaleInvChiSq chisq_;
    Uniform uniform_;
    uvec snp_tracker;
    uvec fold_;
};

template <>
class GeneticSampler<BayesBpi>
    : public BaseGeneticSampler<BayesBpi, GeneticSampler<BayesBpi>>
{
   public:
    GeneticSampler(std::mt19937_64& random_generator, BayesBpi& model);
    void sample_impl(BayesBpi& model, dvec& y_adj);

   private:
    ScaleInvChiSq chisq_;
    Uniform uniform_;
    uvec snp_tracker;
    uvec fold_;
};

template <>
class GeneticSampler<BayesC>
    : public BaseGeneticSampler<BayesC, GeneticSampler<BayesC>>
{
   public:
    GeneticSampler(std::mt19937_64& random_generator, BayesC& model);

    void sample_impl(BayesC& model, dvec& y_adj);

   private:
    ScaleInvChiSq chisq_;
    Uniform uniform_;
    uvec snp_tracker;
    uvec fold_;
};

template <>
class GeneticSampler<BayesCpi>
    : public BaseGeneticSampler<BayesCpi, GeneticSampler<BayesCpi>>
{
   public:
    GeneticSampler(std::mt19937_64& random_generator, BayesCpi& model);

    void sample_impl(BayesCpi& model, dvec& y_adj);

   private:
    ScaleInvChiSq chisq_;
    Uniform uniform_;
    uvec snp_tracker;
    uvec fold_;
};

}  // namespace gelex

#include "gelex/estimator/gibbs/bayes_sample_policy_impl.h"
