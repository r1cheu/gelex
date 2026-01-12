#pragma once

#include <Eigen/Dense>

namespace gelex
{
void constrain(Eigen::Ref<Eigen::VectorXd> varcmp, double y_variance);
}  // namespace gelex
