#include "gelex/predictor/bayes/predictor.h"

#include "gelex/estimator/bayes/result.h"
#include "gelex/model/bayes/model.h"
#include "gelex/utils/utils.h"
namespace gelex
{

BayesPredictor::BayesPredictor(
    const BayesModel& model,
    const MCMCResult& result)
    : formula_(model.formula()), fixed_effects_(result.fixed.mean)
{
    for (size_t i = 0; i < result.random.size(); ++i)
    {
        random_effects_.emplace_back(result.random[i].coeff.mean);
    }
    for (size_t i = 0; i < result.genetic.size(); ++i)
    {
        genetic_effects_.emplace_back(result.snp_eff.col(i));
    }

    for (const auto& genotype : model.genetic())
    {
        genetic_means_.emplace_back(genotype.mean);
        genetic_stds_.emplace_back(genotype.stddev);
    }
};

arma::dvec BayesPredictor::predict(
    const arma::dmat& fixed_design,
    const arma::dcube& random_design,
    arma::dcube& genetic_design) const
{
    arma::dvec prediction = fixed_design * fixed_effects_;

    for (size_t i = 0; i < random_effects_.size(); ++i)
    {
        prediction += random_design.slice(i) * random_effects_[i];
    }

    standardize_genotype(genetic_design);
    for (size_t i = 0; i < genetic_effects_.size(); ++i)
    {
        prediction += genetic_design.slice(i) * genetic_effects_[i];
    }

    return prediction;
}

void BayesPredictor::standardize_genotype(arma::dcube& genotype) const
{
    for (size_t i = 0; i < genotype.n_slices; ++i)
    {
        const auto& genetic_mean = genetic_means_[i];
        const auto& genetic_std = genetic_stds_[i];
#pragma omp parallel for default(none) \
    shared(genotype, i, genetic_mean, genetic_std)
        for (size_t j = 0; j < genotype.n_cols; ++j)
        {
            genotype.col(j)
                = (genotype.col(j) - genetic_mean.at(j)) / genetic_std.at(j);
        }
    };
}

}  // namespace gelex
