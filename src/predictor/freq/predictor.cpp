#include "gelex/predictor/freq/predictor.h"

#include <armadillo>

#include "gelex/data/bed_reader.h"
#include "gelex/model/freq/model.h"

namespace gelex
{
using arma::dmat;
using arma::dvec;
using arma::rowvec;

GBLUPPredictor::GBLUPPredictor(
    std::string_view train_bed,
    const std::vector<std::string>& train_id_order,
    const GBLUP& model)
    : train_bed_(train_bed),
      beta_(model.fixed().beta),
      train_id_order_(train_id_order)
{
    for (const auto& effect : model.random())
    {
        if (effect.type != effect_type::residual)
        {
            blups_.emplace_back(BLUP{effect.level_solutions, effect.name});
        }
    }
}

GBLUPPredictor::GBLUPPredictor(
    std::string_view train_bed,
    const std::vector<std::string>& train_id_order,
    const dvec& beta,
    const std::vector<BLUP>& blups)
    : train_bed_(train_bed),
      beta_(beta),
      blups_(blups),
      train_id_order_(train_id_order)
{
}

void GBLUPPredictor::add_cross_grm(
    std::string_view method,
    dvec p_major,
    double scale_factor,
    size_t chunk_size)
{
    cross_grms_.emplace_back(
        train_bed_,
        std::move(p_major),
        scale_factor,
        chunk_size,
        train_id_order_);
}

dmat GBLUPPredictor::compute_fixed_effects(const dvec& covariates) const
{
    return covariates * beta_;
}

dmat GBLUPPredictor::compute_random_effects(std::string_view test_bed)
{
    auto num_random_effects = cross_grms_.size();

    dmat RandomEffects;

    for (size_t i = 0; i < num_random_effects; ++i)
    {
        dmat new_k = cross_grms_[i].compute(test_bed);
        if (i == 0)
        {
            RandomEffects.zeros(new_k.n_rows, num_random_effects);
        }
        RandomEffects.unsafe_col(i) = new_k * blups_[i].level_solutions;
    }
    return RandomEffects;
}

}  // namespace gelex
