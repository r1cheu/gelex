#include "../src/types/bayes_effects.h"

#include <utility>

#include <Eigen/Core>

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
    std::optional<std::vector<std::string>> levels,
    MatrixXd&& design_matrix)
    : design_matrix(std::move(design_matrix)), levels(std::move(levels))
{
    cols_norm = this->design_matrix.colwise().squaredNorm();
}

FixedState::FixedState(const FixedEffect& effect)
    : coeffs(VectorXd::Zero(effect.design_matrix.cols())) {};

RandomEffect::RandomEffect(
    std::optional<std::vector<std::string>> levels,
    MatrixXd&& design_matrix)
    : design_matrix(std::move(design_matrix)), levels(std::move(levels))
{
    cols_norm = this->design_matrix.colwise().squaredNorm();
}

RandomState::RandomState(const RandomEffect& effect)
    : coeffs(VectorXd::Zero(effect.design_matrix.cols())),
      variance{effect.init_variance}
{
}

// BaseMarkerEffect constructors are inline in header

DominantState::DominantState(const DominantEffect& effect)
    : GeneticState{effect},
      ratios(VectorXd::Zero(get_cols(effect.design_matrix))),
      ratio_mean(effect.ratio_mean),
      ratio_variance(effect.ratio_variance)
{
    if (effect.init_pi)
    {
        tracker = VectorXi::Zero(get_cols(effect.design_matrix));
        pi
            = {effect.init_pi.value(),
               Eigen::VectorXi::Zero(effect.init_pi->size())};
    };
};

}  // namespace bayes
}  // namespace gelex
