#pragma once
#include <cstdint>
#include <memory>
#include <random>

#include <armadillo>

#include "gelex/dist.h"
#include "gelex/estimator/gibbs/base.h"
#include "gelex/estimator/gibbs/bayes_sample_policy.h"
#include "gelex/model/bayes_model.h"

namespace gelex
{
class MCMC
{
   public:
    MCMC(uint64_t n_iter, uint64_t n_burnin, uint64_t n_thin, uint64_t seed);

    template <GeneticPolicy Genetic>
    void Run(BayesianModel<Genetic>& model)
    {
        model.SetPriors();
        Initization(model);
        GeneticSampler genetic_sampler(gen_, model);

        for (uint64_t i = 0; i < n_iter_; ++i)
        {
            SampleMu(model);
            SampleFixedEffects(model);
            SampleEnvironmentalEffect(model);

            genetic_sampler.Sample(model, y_adj_);
            SampleErrorVariance(model);

            if (i >= n_burnin_ && i % n_thin_ == 0)
            {
                StoreSample(model);
            }
        }
    }

   private:
    auto ComputeNumRecords() const
    {
        return static_cast<arma::uword>((n_iter_ - n_burnin_) / n_thin_);
    }

    template <GeneticPolicy Genetic>
    void SampleMu(BayesianModel<Genetic>& model)
    {
        double mu = model.mu();                             // NOLINT
        const auto n = static_cast<double>(y_adj_.n_elem);  // NOLINT
        double mu_adj
            = -(*normal_)(arma::sum(y_adj_) / n, sqrt(model.sigma_e() / n));
        mu -= mu_adj;
        daxpy_auto(y_adj_, ones_, mu_adj);
        model.set_mu(mu);
    }

    template <GeneticPolicy Genetic>
    void SampleFixedEffects(BayesianModel<Genetic>& model)
    {
        if (!model.HasBeta())
        {
            return;
        }

        SampleEffect(
            *normal_,
            model.beta(),
            y_adj_,
            *(model.design_mat_beta()),
            model.beta_col_norm2(),
            model.sigma_e());
    }

    template <GeneticPolicy Genetic>
    void SampleEnvironmentalEffect(BayesianModel<Genetic>& model)
    {
        if (!model.HasEnv())
        {
            return;
        }

        SampleEffect(
            *normal_,
            model.r(),
            y_adj_,
            *(model.design_mat_r()),
            model.r_col_norm2(),
            model.sigma_e(),
            model.sigma_r());

        const double ssq = arma::dot(model.r(), model.r());
        model.set_sigma_r((*sigma_r_sampler_)(ssq));
    }

    template <GeneticPolicy Genetic>
    void SampleErrorVariance(BayesianModel<Genetic>& model)
    {
        const double ssq = arma::dot(y_adj_, y_adj_);
        model.set_sigma_e((*sigma_e_sampler_)(ssq));
    }

    template <GeneticPolicy Genetic>
    void Initization(BayesianModel<Genetic>& model)
    {
        store_idx_ = 0;
        y_adj_ = model.phenotype() - arma::mean(model.phenotype());
        ones_ = arma::ones<dvec>(y_adj_.n_elem);
        n_records_ = ComputeNumRecords();

        normal_ = std::make_unique<Normal>(gen_, 0.0, 1.0);
        if (model.HasEnv())
        {
            sigma_r_sampler_ = std::make_unique<ScaleInvChiSq>(
                gen_,
                model.priors().sigma_r().nu,
                static_cast<double>(model.r().n_cols),
                model.priors().sigma_r().s2);
        }

        sigma_e_sampler_ = std::make_unique<ScaleInvChiSq>(
            gen_,
            model.priors().sigma_e().nu,
            static_cast<double>(y_adj_.n_elem),
            model.priors().sigma_e().s2);
    }

    template <GeneticPolicy Genetic>
    void InitStorage(BayesianModel<Genetic>& model)
    {
        mu_store = arma::zeros<dvec>(n_records_);

        if (model.HasBeta())
        {
            beta_store = arma::zeros<dmat>(model.beta().n_elem, n_records_);
        }

        a_store = arma::zeros<dmat>(model.a().n_elem, n_records_);

        if (model.HasEnv())
        {
            r_store = arma::zeros<dmat>(model.r().n_elem, n_records_);
            sigma_r_store = arma::zeros<dvec>(n_records_);
        }

        if constexpr (std::is_same_v<
                          typename BayesianModel<Genetic>::sigma_t,
                          double>)
        {
            sigma_a_store = arma::zeros<dmat>(1, n_records_);
        }
        else
        {
            sigma_a_store
                = arma::zeros<dmat>(model.sigma_a().n_elem, n_records_);
        }

        sigma_e_store = arma::zeros<dvec>(n_records_);
    }

    template <GeneticPolicy Genetic>
    void StoreSample(BayesianModel<Genetic>& model)
    {
        mu_store.at(store_idx_) = model.mu();

        if (model.HasBeta())
        {
            beta_store.col(store_idx_) = model.beta();
        }

        a_store.col(store_idx_) = model.a();

        if (model.HasEnv())
        {
            r_store.col(store_idx_) = model.r();
            sigma_r_store.at(store_idx_) = model.sigma_r();
        }

        sigma_a_store.col(store_idx_) = model.sigma_a();
        sigma_e_store.at(store_idx_) = model.sigma_e();

        ++store_idx_;
    };

    // --- random number generator ---
    uint64_t n_iter_;
    uint64_t n_burnin_;
    uint64_t n_thin_;
    uint64_t n_records_;
    uint64_t store_idx_{};

    std::unique_ptr<Normal> normal_;
    std::unique_ptr<ScaleInvChiSq> sigma_e_sampler_;
    std::unique_ptr<ScaleInvChiSq> sigma_r_sampler_;
    std::mt19937_64 gen_;

    dvec y_adj_;
    dvec ones_;

    dvec mu_store;
    dmat beta_store;
    dmat a_store;
    dmat r_store;

    dmat sigma_a_store;
    dvec sigma_r_store;
    dvec sigma_e_store;
};

}  // namespace gelex
