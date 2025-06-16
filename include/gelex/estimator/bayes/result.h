#pragma once
#include <string>
#include <vector>

#include <armadillo>
#include "gelex/estimator/bayes/samples.h"

namespace gelex
{

class MCMCSamples;
class BayesModel;

struct PosteriorStats
{
    arma::dvec means;
    arma::dvec stds;
    arma::dvec medians;
    arma::dvec q5s;
    arma::dvec q95s;
    arma::dvec n_effs;
    arma::dvec r_hats;

    explicit operator bool() const { return means.n_elem > 0; }

    double mean(size_t index) const { return means.at(index); }
    double std(size_t index) const { return stds.at(index); }
    double median(size_t index) const { return medians.at(index); }
    double q5(size_t index) const { return q5s.at(index); }
    double q95(size_t index) const { return q95s.at(index); }
    double n_eff(size_t index) const { return n_effs.at(index); }
    double r_hat(size_t index) const { return r_hats.at(index); }
};

struct PosteriorGroup
{
    std::vector<PosteriorStats> coeffs;
    std::vector<PosteriorStats> sigmas;
    std::vector<std::string> names;
};

struct MCMCResult
{
    PosteriorStats mu;
    PosteriorStats fixed;
    PosteriorGroup random;
    PosteriorGroup genetic;
    PosteriorStats residual;
    PosteriorStats h2;
};

MCMCResult compute_mcmc_result(
    const MCMCSamples& samples,
    const BayesModel& model);

PosteriorStats compute_effect_result(const arma::dvec& samples);
PosteriorStats compute_effect_result(const arma::dmat& samples);
//
// template <typename EffectManager>
// void process_posterior(
//     const SampleGroup& samples,
//     const EffectManager& effect,
//     PosteriorGroup& result)
// {
//     auto n_eff = effect.size();
//     result.coeffs.reserve(n_eff);
//     result.sigmas.reserve(n_eff);
//     result.names.reserve(n_eff);
//
//     for (size_t i = 0; i < n_eff; ++i)
//     {
//         result.coeffs.push_back(compute_effect_result(samples.coeffs[i]));
//         result.sigmas.push_back(compute_effect_result(samples.sigmas[i]));
//         result.names.push_back(effect[i].name);
//     }
// }

}  // namespace gelex
