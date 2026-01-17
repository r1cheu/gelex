#ifndef GELEX_OPTIM_CONVERGENCE_CHECKER_H_
#define GELEX_OPTIM_CONVERGENCE_CHECKER_H_

#include <Eigen/Dense>

namespace gelex
{
class ConvergenceChecker
{
   public:
    explicit ConvergenceChecker(double tol = 1e-8) : tol_{tol} {}

    bool is_converged(
        const Eigen::Ref<const Eigen::VectorXd>& new_sigma,
        double new_loglike);

    void clear()
    {
        converged_ = false;
        old_sigma_.resize(0);
        old_loglike_ = 0.0;
    }

   private:
    double tol_;
    bool converged_{false};
    Eigen::VectorXd old_sigma_;
    double old_loglike_{};

    void update(const Eigen::Ref<const Eigen::VectorXd>& sigma, double loglike)
    {
        old_sigma_ = sigma;
        old_loglike_ = loglike;
    }

    double compute_sigma_diff(
        const Eigen::Ref<const Eigen::VectorXd>& new_sigma) const;
    double compute_loglike_diff(double new_loglike) const;
};
}  // namespace gelex

#endif  // GELEX_OPTIM_CONVERGENCE_CHECKER_H_
