#include "bayes_effects.h"

#include <fstream>
#include <unordered_set>
#include <utility>

#include "Eigen/Core"

#include "../src/data/loader.h"
#include "gelex/model/effects.h"

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
    std::vector<std::string>&& levels,
    MatrixXd&& design_matrix)
    : names(std::move(names)),
      levels(std::move(levels)),
      design_matrix(std::move(design_matrix))
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
      sigma(effect.sigma)  // Use dvec for consistency with other effects
{
}

AdditiveEffect::AdditiveEffect(
    std::string&& name,
    const std::string& genotype_bin,
    const std::string& bim_file_path,
    BayesAlphabet type,
    VectorXd&& sigma,
    VectorXd&& pi)
    : name(std::move(name)),
      pi(std::move(pi)),
      sigma(std::move(sigma)),
      type(type),
      design_matrix(genotype_bin),
      bim_file_path(bim_file_path),
      mono_indices(detail::read_mono_indices(genotype_bin))
{
    cols_norm = this->design_matrix.mat.colwise().squaredNorm();
}

AdditiveStatus::AdditiveStatus(const AdditiveEffect& effect)
    : u(VectorXd::Zero(effect.design_matrix.mat.rows())),
      coeff(VectorXd::Zero(effect.design_matrix.mat.cols())),
      tracker(VectorXi::Zero(effect.design_matrix.mat.cols())),
      type(effect.type),
      pi{effect.pi, Eigen::VectorXi::Zero(effect.pi.size())},
      sigma(effect.sigma) {};

DominantEffect::DominantEffect(
    std::string&& name,
    const std::string& genotype_bin,
    double prior_mean,
    double prior_var)
    : name(std::move(name)),
      design_matrix(genotype_bin),
      prior_mean(prior_mean),
      prior_var(prior_var),
      mono_indices(detail::read_mono_indices(genotype_bin))
{
    cols_norm = this->design_matrix.mat.colwise().squaredNorm();
}

DominantStatus::DominantStatus(const DominantEffect& effect)
    : coeff(VectorXd::Zero(effect.design_matrix.mat.cols())),
      u(VectorXd::Zero(effect.design_matrix.mat.cols()))
{
}

bool AdditiveEffect::is_monomorphic(Eigen::Index snp_index) const
{
    return mono_indices.contains(snp_index);
}

Index AdditiveEffect::num_mono() const
{
    return static_cast<Index>(mono_indices.size());
}

bool DominantEffect::is_monomorphic(Eigen::Index snp_index) const
{
    return mono_indices.contains(snp_index);
}

Index DominantEffect::num_mono() const
{
    return static_cast<Index>(mono_indices.size());
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
