#include "gelex/model/bayes/effects.h"
#include "gelex/model/effects.h"
#include "gelex/utils/utils.h"

namespace gelex
{
using arma::dmat;
using arma::dvec;
using arma::uvec;

namespace bayes
{

FixedEffect::FixedEffect(
    std::vector<std::string>&& names,
    std::vector<std::string>&& levels,
    dmat&& design_matrix)
    : names(std::move(names)),
      levels(std::move(levels)),
      design_matrix(std::move(design_matrix))
{
    cols_norm = sum_square(this->design_matrix);
}

FixedEffectState::FixedEffectState(size_t n_coeff)
    : coeff(n_coeff, arma::fill::zeros) {};

RandomEffect::RandomEffect(std::string&& name, dmat&& design_matrix)
    : name(std::move(name)), design_matrix(std::move(design_matrix))
{
    cols_norm = sum_square(this->design_matrix);
};

RandomEffectState::RandomEffectState(size_t n_coeff, const dvec& init_sigma)
    : coeff(n_coeff, arma::fill::zeros), sigma{init_sigma}
{
}

GeneticEffect::GeneticEffect(
    std::string&& name,
    dmat&& design_matrix,
    BayesAlphabet type,
    dvec&& sigma,
    dvec&& pi)
    : name(std::move(name)),
      design_matrix(std::move(design_matrix)),
      pi(std::move(pi)),
      sigma(std::move(sigma)),
      type(type)
{
    std::tie(mean, stddev) = standradize(this->design_matrix);
}

GeneticEffectState::GeneticEffectState(
    BayesAlphabet type,
    size_t n_individual,
    size_t n_coeff,
    const dvec& pi_prop,
    const dvec& sigma)
    : u(n_individual, arma::fill::zeros),
      coeff(n_coeff, arma::fill::zeros),
      type(type),
      pi{pi_prop, uvec(pi_prop.n_elem, arma::fill::zeros)},
      sigma(sigma) {};
std::vector<GeneticEffectState> create_thread_states(
    const GeneticEffectManager& designs)
{
    std::vector<GeneticEffectState> states;
    states.reserve(designs.size());
    for (const auto& design : designs.effects())
    {
        states.emplace_back(
            design.type,
            design.design_matrix.n_rows,
            design.design_matrix.n_cols,
            design.pi,
            design.sigma);
    }
    return states;
}

std::vector<RandomEffectState> create_thread_states(
    const RandomEffectManager& designs)
{
    std::vector<RandomEffectState> states;
    states.reserve(designs.size());
    for (const auto& design : designs.effects())
    {
        states.emplace_back(design.design_matrix.n_cols, design.sigma);
    }
    return states;
}
}  // namespace bayes

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
