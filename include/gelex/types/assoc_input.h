#ifndef GELEX_TYPES_ASSOC_INPUT_H_
#define GELEX_TYPES_ASSOC_INPUT_H_

#include <Eigen/Core>

namespace gelex
{

struct AssocInput
{
    AssocInput(
        Eigen::Index chunk_size,
        Eigen::MatrixXd V_inv_in,
        Eigen::VectorXd y_adj)
        : Z(Eigen::MatrixXd::Zero(V_inv_in.rows(), chunk_size)),
          V_inv(std::move(V_inv_in)),
          V_inv_y(std::move(y_adj)),
          W(Eigen::MatrixXd::Zero(V_inv.rows(), chunk_size))
    {
    }
    Eigen::MatrixXd Z;  // SNP matrix
    Eigen::MatrixXd V_inv;
    Eigen::VectorXd V_inv_y;
    Eigen::MatrixXd W;  // Intermediate buffer for V^{-1} Z
};

struct AssocOutput
{
    explicit AssocOutput(Eigen::Index chunk_size)
        : beta(Eigen::VectorXd::Zero(chunk_size)),
          se(Eigen::VectorXd::Zero(chunk_size)),
          stats(Eigen::VectorXd::Zero(chunk_size)),
          p_value(Eigen::VectorXd::Zero(chunk_size)),
          zt_v_inv_r(Eigen::VectorXd::Zero(chunk_size)),
          zt_v_inv_z(Eigen::VectorXd::Zero(chunk_size))
    {
    }
    Eigen::VectorXd beta;
    Eigen::VectorXd se;
    Eigen::VectorXd stats;
    Eigen::VectorXd p_value;

    Eigen::VectorXd zt_v_inv_r;  // Z^t V^{-1} residual
    Eigen::VectorXd zt_v_inv_z;  // Z^t V^{-1} Z
};

}  // namespace gelex

#endif  // GELEX_TYPES_ASSOC_INPUT_H_
