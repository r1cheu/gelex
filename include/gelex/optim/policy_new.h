#ifndef GELEX_OPTIM_POLICY_NEW_H_
#define GELEX_OPTIM_POLICY_NEW_H_

#include <Eigen/Core>

namespace gelex
{

class FreqModel;
class FreqState;
class OptimizerState;

struct EMPolicy
{
    static auto apply(
        const FreqModel& model,
        const FreqState& state,
        OptimizerState& opt_state) -> Eigen::VectorXd;
};

struct AIPolicy
{
    static auto apply(
        const FreqModel& model,
        const FreqState& state,
        OptimizerState& opt_state) -> Eigen::VectorXd;
};

}  // namespace gelex

#endif  // GELEX_OPTIM_POLICY_NEW_H_
