#include "../src/types/freq_effect.h"

namespace gelex::freq
{

FixedState::FixedState(const gelex::FixedEffect& effect)
    : coeff(Eigen::VectorXd::Zero(effect.X.cols())),
      se(Eigen::VectorXd::Zero(effect.X.cols()))
{
}

RandomState::RandomState(const RandomEffect& effect)
    : name(effect.name),
      blup(Eigen::VectorXd::Zero(effect.levels.size())),
      variance(0.0)
{
}

GeneticState::GeneticState(const GeneticEffect& effect)
    : name(effect.name),
      ebv(Eigen::VectorXd::Zero(effect.K.rows())),
      variance(0.0)
{
}

}  // namespace gelex::freq
