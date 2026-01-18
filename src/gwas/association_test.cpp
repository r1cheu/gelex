#include "gelex/gwas/association_test.h"
#include <mkl_vml_functions.h>

#include <Eigen/Dense>

namespace gelex::gwas
{

namespace
{

void chi2_test(
    Eigen::Ref<Eigen::VectorXd> wald_stats,
    Eigen::Ref<Eigen::VectorXd> p_values)
{
    int n = static_cast<int>(wald_stats.size());
    vdSqrt(n, wald_stats.data(), p_values.data());
    cblas_dscal(n, -1.0, p_values.data(), 1);
    vdCdfNorm(n, p_values.data(), p_values.data());
    cblas_dscal(n, 2.0, p_values.data(), 1);
}
}  // namespace

void wald_test(const AssocInput& input, AssocOutput& output)
{
    output.zt_v_inv_r = (input.Z.transpose() * input.V_inv_y);
    output.zt_v_inv_z = (input.Z.transpose() * input.V_inv * input.Z)
                            .diagonal();  // should be optimized by Eigen

    output.beta = (output.zt_v_inv_r.array() / output.zt_v_inv_z.array());
    output.se = (1.0 / output.zt_v_inv_z.array()).sqrt();
    output.stats = (output.beta.array() / output.se.array()).square();
    chi2_test(output.stats, output.p_value);
}
}  // namespace gelex::gwas
