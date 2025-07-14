#include "gelex/predictor/freq/predictor.h"

#include <armadillo>

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
      beta_(model.fixed().coeff),
      train_id_order_(train_id_order),
      proj_y_(model.proj_y_)
{
    for (const auto& effect : model.random_)
    {
        blups_.emplace_back(effect.name, effect.design_matrix, effect.sigma);
    }
    for (const auto& effect : model.genetic_)
    {
        blups_.emplace_back(effect.name, effect.design_matrix, effect.sigma);
    }
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

dvec GBLUPPredictor::compute_random_effects(
    std::string_view name,
    const dmat& design_matrix) const
{
    auto it = std::find_if(
        blups_.begin(),
        blups_.end(),
        [&name](const Effect& blup) { return blup.name == name; });
    if (it != blups_.end())
    {
        return design_matrix * it->design_matrix.t() * proj_y_ * it->sigma;
    }

    throw std::runtime_error("Random effect not found: " + std::string(name));
}

dmat GBLUPPredictor::compute_genetic_effects(std::string_view test_bed)
{
    auto num_effects = cross_grms_.size();

    dmat genetic_effects;

    for (size_t i = 0; i < num_effects; ++i)
    {
        dmat new_k = cross_grms_[i].compute(test_bed);
        if (i == 0)
        {
            genetic_effects.zeros(new_k.n_rows, num_effects);
        }
        genetic_effects.unsafe_col(i)
            = new_k * blups_[i].design_matrix.t() * proj_y_ * blups_[i].sigma;
    }
    return genetic_effects;
}

}  // namespace gelex
