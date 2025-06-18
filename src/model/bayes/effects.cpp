#include "gelex/model/bayes/effects.h"
#include "gelex/model/effects.h"

namespace gelex
{
using arma::dmat;
using arma::dvec;
using arma::uvec;

BaseEffectDesign::BaseEffectDesign(dmat&& design_mat_)
    : design_mat(std::move(design_mat_))
{
    cols_norm = sum_square(design_mat);
};

BaseEffectState::BaseEffectState(size_t n_coeff)
    : coeff(n_coeff, arma::fill::zeros)
{
}

FixedEffectDesign::FixedEffectDesign(
    std::vector<std::string>&& names_,
    std::vector<std::string>&& levels_,
    dmat&& design_mat_)
    : names(std::move(names_)),
      levels(std::move(levels_)),
      BaseEffectDesign(std::move(design_mat_))
{
}

RandomEffectDesign::RandomEffectDesign(std::string&& name_, dmat&& design_mat_)
    : name(std::move(name_)), BaseEffectDesign(std::move(design_mat_)) {};

RandomEffectState::RandomEffectState(size_t n_coeff, const dvec& init_sigma)
    : BaseEffectState(n_coeff), sigma{init_sigma}
{
}

GeneticEffectDesign::GeneticEffectDesign(
    std::string&& name_,
    dmat&& design_mat_,
    BayesAlphabet type_,
    dvec&& sigma_,
    dvec&& pi_)
    : RandomEffectDesign(std::move(name_), std::move(design_mat_)),
      pi(std::move(pi_)),
      sigma(std::move(sigma_)),
      type(type_)
{
}

GeneticEffectState::GeneticEffectState(
    size_t n_individual,
    size_t n_coeff,
    const dvec& pi_prop,
    const dvec& sigma_)
    : u(n_individual, arma::fill::zeros),
      BaseEffectState(n_coeff),
      pi{pi_prop, uvec(pi_prop.n_elem, arma::fill::zeros)},
      sigma(sigma_) {};

std::vector<GeneticEffectState> create_thread_states(
    const GeneticEffectDesignManager& designs)
{
    std::vector<GeneticEffectState> states;
    states.reserve(designs.size());
    for (const auto& design : designs.effects())
    {
        states.emplace_back(
            design.design_mat.n_rows,
            design.design_mat.n_cols,
            design.pi,
            design.sigma);
    }
    return states;
}

std::vector<RandomEffectState> create_thread_states(
    const RandomEffectDesignManager& designs)
{
    std::vector<RandomEffectState> states;
    states.reserve(designs.size());
    for (const auto& design : designs.effects())
    {
        states.emplace_back(design.design_mat.n_cols, design.sigma);
    }
    return states;
}

dvec compute_cols_var(const dmat& mat)
{
    dvec out(mat.n_cols);
#pragma omp parallel for default(none) shared(mat, out)
    for (size_t i = 0; i < mat.n_cols; ++i)
    {
        out.at(i) = arma::var(mat.unsafe_col(i));
    }
    return out;
}

}  // namespace gelex
