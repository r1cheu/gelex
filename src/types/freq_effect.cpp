#include "../src/types/freq_effect.h"

#include <Eigen/Core>

namespace gelex::freq
{

FixedState::FixedState(const gelex::FixedEffect& effect)
    : coeff(Eigen::VectorXd::Zero(effect.X.cols())),
      se(Eigen::VectorXd::Zero(effect.X.cols()))
{
}

RandomState::RandomState(const RandomEffect& effect)
    : name(effect.name),
      blup(
          Eigen::VectorXd::Zero(
              static_cast<Eigen::Index>(effect.levels.size())))
{
}

GeneticState::GeneticState(const GeneticEffect& effect)
    : name(effect.name), ebv(Eigen::VectorXd::Zero(effect.K.rows()))
{
}

}  // namespace gelex::freq
