#ifndef GELEX_TYPES_ASSOC_INPUT_H_
#define GELEX_TYPES_ASSOC_INPUT_H_

#include <Eigen/Core>

namespace gelex
{

struct AssocInput
{
    Eigen::MatrixXd V_inv;
    Eigen::VectorXd y_adj;
};

}  // namespace gelex

#endif  // GELEX_TYPES_ASSOC_INPUT_H_
