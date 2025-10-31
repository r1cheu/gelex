#pragma once

#include <armadillo>
#include "gelex/model/bayes/model.h"
#include "gelex/types/mcmc_results.h"
namespace gelex
{

class BayesPredictor
{
   public:
    BayesPredictor(const BayesModel& model, const MCMCResult& result);

    arma::dvec predict(
        const arma::dmat& fixed_design,
        const arma::dcube& random_design,
        arma::dcube& genetic_design) const;

   private:
    void standardize_genotype(arma::dcube& genotype) const;

    std::string formula_;
    arma::dvec fixed_effects_;
    std::vector<arma::dvec> random_effects_;
    std::vector<arma::dvec> genetic_effects_;

    std::vector<arma::dvec> genetic_means_;
    std::vector<arma::dvec> genetic_stds_;
};
}  // namespace gelex
