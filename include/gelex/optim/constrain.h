#ifndef GELEX_OPTIM_CONSTRAIN_H_
#define GELEX_OPTIM_CONSTRAIN_H_

#include <Eigen/Dense>

namespace gelex
{
void constrain(Eigen::Ref<Eigen::VectorXd> varcmp, double y_variance);
}  // namespace gelex

#endif  // GELEX_OPTIM_CONSTRAIN_H_
