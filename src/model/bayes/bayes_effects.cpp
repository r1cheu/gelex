#include "bayes_effects.h"

#include <utility>

#include "Eigen/Core"

namespace gelex
{
using Eigen::Ref;

using Eigen::Index;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

namespace bayes
{

FixedEffect::FixedEffect(
    std::vector<std::string>&& names,
    MatrixXd&& design_matrix)
    : names(std::move(names)), design_matrix(std::move(design_matrix))
{
    cols_norm = this->design_matrix.colwise().squaredNorm();
}

FixedStatus::FixedStatus(Index n_coeff) : coeff(VectorXd::Zero(n_coeff)) {};

RandomEffect::RandomEffect(std::string&& name, MatrixXd&& design_matrix)
    : name(std::move(name)), design_matrix(std::move(design_matrix))
{
    cols_norm = this->design_matrix.colwise().squaredNorm();
};

RandomStatus::RandomStatus(const RandomEffect& effect)
    : coeff(VectorXd::Zero(effect.design_matrix.cols())),
      effect_variance(
          effect
              .effect_variance)  // Use dvec for consistency with other effects
{
}

AdditiveEffect::AdditiveEffect(
    std::string&& name,
    GenotypeMap&& design_matrix,
    VectorXd&& sigma,
    VectorXd&& pi)
    : name(std::move(name)),
      pi(std::move(pi)),
      marker_variance(std::move(sigma)),
      design_matrix(std::move(design_matrix))
{
    cols_norm = this->design_matrix.matrix().colwise().squaredNorm();
}

AdditiveStatus::AdditiveStatus(const AdditiveEffect& effect)
    : u(VectorXd::Zero(effect.design_matrix.rows())),
      coeff(VectorXd::Zero(effect.design_matrix.cols())),
      tracker(VectorXi::Zero(effect.design_matrix.cols())),
      pi{effect.pi, Eigen::VectorXi::Zero(effect.pi.size())},
      marker_variance(effect.marker_variance) {};

DominantEffect::DominantEffect(
    std::string&& name,
    GenotypeMap&& design_matrix,
    double prior_mean,
    double prior_var)
    : name(std::move(name)),
      design_matrix(std::move(design_matrix)),
      prior_mean(prior_mean),
      prior_var(prior_var)
{
    cols_norm = this->design_matrix.matrix().colwise().squaredNorm();
}

DominantStatus::DominantStatus(const DominantEffect& effect)
    : coeff(VectorXd::Zero(effect.design_matrix.cols())),
      u(VectorXd::Zero(effect.design_matrix.rows()))
{
}

bool AdditiveEffect::is_monomorphic(Eigen::Index snp_index) const
{
    return design_matrix.is_monomorphic(snp_index);
}

Index AdditiveEffect::num_mono() const
{
    return design_matrix.num_mono();
}

bool DominantEffect::is_monomorphic(Eigen::Index snp_index) const
{
    return design_matrix.is_monomorphic(snp_index);
}

Index DominantEffect::num_mono() const
{
    return design_matrix.num_mono();
}

std::vector<RandomStatus> create_chain_states(
    const RandomEffectManager& effects)
{
    std::vector<RandomStatus> states;
    states.reserve(effects.size());
    for (const auto& effect : effects.effects())
    {
        states.emplace_back(effect);
    }
    return states;
}
}  // namespace bayes
}  // namespace gelex
