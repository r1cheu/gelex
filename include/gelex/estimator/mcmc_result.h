#pragma once
#include <armadillo>
#include <string>
#include <vector>

namespace gelex
{

class MCMCStorage;
class Bayes;

struct ParameterResult
{
    double mean;
    double std;
    double median;
    double q5;     // 5th percentile
    double q95;    // 95th percentile
    double n_eff;  // effective sample size
    double r_hat;  // potential scale reduction factor
};

struct EffectResult
{
    std::vector<ParameterResult> parameters;

    double mean(size_t index) const { return parameters.at(index).mean; }
    double std(size_t index) const { return parameters.at(index).std; }
    double median(size_t index) const { return parameters.at(index).median; }
    double q5(size_t index) const { return parameters.at(index).q5; }
    double q95(size_t index) const { return parameters.at(index).q95; }
    double n_eff(size_t index) const { return parameters.at(index).n_eff; }
    double r_hat(size_t index) const { return parameters.at(index).r_hat; }
    arma::dvec means() const { return extract_values(&ParameterResult::mean); }
    arma::dvec stds() const { return extract_values(&ParameterResult::std); }
    arma::dvec medians() const
    {
        return extract_values(&ParameterResult::median);
    }
    arma::dvec q5s() const { return extract_values(&ParameterResult::q5); }
    arma::dvec q95s() const { return extract_values(&ParameterResult::q95); }
    arma::dvec n_effs() const
    {
        return extract_values(&ParameterResult::n_eff);
    }
    arma::dvec r_hats() const
    {
        return extract_values(&ParameterResult::r_hat);
    }

   private:
    template <typename MemberPtr>
    arma::dvec extract_values(MemberPtr member) const
    {
        arma::dvec result(parameters.size());
        for (size_t i = 0; i < parameters.size(); ++i)
        {
            result.at(i) = parameters[i].*member;
        }
        return result;
    }
};

struct MCMCResult
{
    EffectResult mu;
    EffectResult fixed;
    std::vector<EffectResult> random;
    std::vector<EffectResult> genetic;
    EffectResult residual;

    std::vector<EffectResult> random_sigma;
    std::vector<EffectResult> genetic_sigma;

    std::vector<std::string> random_names;
    std::vector<std::string> genetic_names;
};

MCMCResult compute_mcmc_result(const MCMCStorage& storage, const Bayes& model);

}  // namespace gelex
