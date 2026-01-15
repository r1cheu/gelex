#ifndef GELEX_OPTIM_OPTIMIZER_STATE_H_
#define GELEX_OPTIM_OPTIMIZER_STATE_H_

#include <Eigen/Core>

namespace gelex
{

class FreqModel;

class OptimizerState
{
   public:
    explicit OptimizerState(const FreqModel& model);

    auto phenotype_var() const -> double { return phenotype_var_; }
    auto num_individuals() const -> Eigen::Index { return num_individuals_; }

    // computed matrices, public for Policy access
    Eigen::MatrixXd v;
    Eigen::MatrixXd proj;
    Eigen::VectorXd proj_y;
    Eigen::MatrixXd tx_vinv_x;
    double logdet_v{};

    // for AI policy
    Eigen::MatrixXd hess_inv;

   private:
    Eigen::Index num_individuals_{};
    double phenotype_var_{};
};

}  // namespace gelex

#endif  // GELEX_OPTIM_OPTIMIZER_STATE_H_
