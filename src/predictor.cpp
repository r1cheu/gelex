#include "gelex/predictor.h"

#include <memory>

#include <armadillo>

#include "gelex/data/cross_grm.h"
#include "gelex/model/linear_mixed_model.h"

namespace gelex
{
using arma::rowvec;
void Predictor::set_cross_grm(
    std::string_view method,
    rowvec center,
    double scale_factor,
    uint64_t chunk_size)
{
    if (method == "add")
    {
        cross_grms_.emplace_back(
            std::make_unique<AddCrossGrm>(
                train_bed_,
                std::move(center),
                scale_factor,
                chunk_size,
                params_.dropped_individuals()));
    }
    else if (method == "dom")
    {
        cross_grms_.emplace_back(
            std::make_unique<DomCrossGrm>(
                train_bed_,
                std::move(center),
                scale_factor,
                chunk_size,
                params_.dropped_individuals()));
    }
    else
    {
        throw std::runtime_error("Unknown method: " + std::string(method));
    }
}

dmat Predictor::ComputeFixedEffects(const dvec& covariates) const noexcept
{
    return covariates * params_.beta();
}

dmat Predictor::ComputeRandomEffects(std::string_view test_bed)
{
    const dvec& proj_y = params_.proj_y();
    const dvec& sigma = params_.sigma();
    auto num_random_effects = cross_grms_.size();

    dmat RandomEffects;

    for (size_t i = 0; i < num_random_effects; ++i)
    {
        dmat new_k = cross_grms_[i]->Compute(test_bed);
        if (i == 0)
        {
            RandomEffects.zeros(new_k.n_rows, num_random_effects);
        }
        RandomEffects.col(i) = new_k * proj_y * sigma[i];
    }
    return RandomEffects;
}

}  // namespace gelex
