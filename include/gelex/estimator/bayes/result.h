#pragma once
#include <armadillo>
#include <vector>
#include "gelex/estimator/bayes/samples.h"

namespace gelex
{
struct PosteriorSummary
{
    explicit PosteriorSummary(size_t n_params)
        : mean(n_params),
          median(n_params),
          stddev(n_params),
          hpdi_low(n_params),
          hpdi_high(n_params),
          ess(n_params),
          rhat(n_params)
    {
    }

    arma::dvec mean;
    arma::dvec median;
    arma::dvec stddev;
    arma::dvec hpdi_low;
    arma::dvec hpdi_high;
    arma::dvec ess;   // Effective sample size
    arma::dvec rhat;  // Gelman-Rubin diagnostic

    std::vector<double> operator[](size_t idx) const
    {
        return std::vector<double>{
            mean.at(idx),
            median.at(idx),
            stddev.at(idx),
            hpdi_low.at(idx),
            hpdi_high.at(idx),
            ess.at(idx),
            rhat.at(idx)};
    }
};

struct PostieriorRandomSummary
{
    explicit PostieriorRandomSummary(size_t n_coeff, size_t n_sigma)
        : coeff(n_coeff), sigma(n_sigma)
    {
    }

    PosteriorSummary coeff;
    PosteriorSummary sigma;
};

struct PostieriorGeneticSummary : PostieriorRandomSummary
{
    PostieriorGeneticSummary(size_t n_coeff, size_t n_sigma, size_t n_pi)
        : PostieriorRandomSummary(n_coeff, n_sigma), pi(n_pi)
    {
    }

    PosteriorSummary heritability{1};
    PosteriorSummary genetic_var{1};
    PosteriorSummary pi;
};

class MCMCResult
{
   public:
    explicit MCMCResult(const MCMCSamples& samples);

    void compute_summary_statistics(
        const MCMCSamples& samples,
        double prob = 0.9);

    PosteriorSummary fixed;
    std::vector<PostieriorRandomSummary> random;
    std::vector<PostieriorGeneticSummary> genetic;
    PosteriorSummary residual;

   private:
    static void compute_summary_statistics(
        PosteriorSummary& summary,
        const arma::dcube& samples,
        double prob = 0.9);
};

}  // namespace gelex
