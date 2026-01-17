#include "gelex/gwas/association_test.h"

#include <cmath>

#include <Eigen/Dense>

namespace gelex::gwas
{

namespace
{

// Lower incomplete gamma function using series expansion
// gamma(a, x) = integral from 0 to x of t^(a-1) * e^(-t) dt
auto lower_incomplete_gamma(double a, double x) -> double
{
    if (x < 0.0)
    {
        return 0.0;
    }
    if (x == 0.0)
    {
        return 0.0;
    }

    // Series expansion: gamma(a,x) = x^a * e^(-x) * sum_{n=0}^inf x^n /
    // (a+1)...(a+n)
    const int max_iter = 200;
    const double eps = 1e-15;

    double sum = 1.0 / a;
    double term = 1.0 / a;

    for (int n = 1; n < max_iter; ++n)
    {
        term *= x / (a + n);
        sum += term;
        if (std::abs(term) < eps * std::abs(sum))
        {
            break;
        }
    }

    return std::pow(x, a) * std::exp(-x) * sum;
}

// Regularized lower incomplete gamma function P(a, x) = gamma(a, x) / Gamma(a)
auto regularized_gamma_p(double a, double x) -> double
{
    if (x < 0.0 || a <= 0.0)
    {
        return 0.0;
    }
    if (x == 0.0)
    {
        return 0.0;
    }

    // For large x, use continued fraction (complementary)
    if (x > a + 1.0)
    {
        // Use Q(a,x) = 1 - P(a,x) via continued fraction
        // Then return 1 - Q
        const int max_iter = 200;
        const double eps = 1e-15;

        double b = x + 1.0 - a;
        double c = 1.0 / 1e-30;
        double d = 1.0 / b;
        double h = d;

        for (int n = 1; n < max_iter; ++n)
        {
            const double an = -n * (n - a);
            b += 2.0;
            d = an * d + b;
            if (std::abs(d) < 1e-30)
            {
                d = 1e-30;
            }
            c = b + an / c;
            if (std::abs(c) < 1e-30)
            {
                c = 1e-30;
            }
            d = 1.0 / d;
            const double del = d * c;
            h *= del;
            if (std::abs(del - 1.0) < eps)
            {
                break;
            }
        }

        const double q = std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
        return 1.0 - q;
    }

    // Use series expansion for small x
    return lower_incomplete_gamma(a, x) / std::tgamma(a);
}

// Chi-squared CDF: P(X <= x) where X ~ chi^2(df)
// chi^2 CDF = P(df/2, x/2) where P is regularized lower incomplete gamma
auto chi_squared_cdf(double x, int df) -> double
{
    if (x <= 0.0 || df <= 0)
    {
        return 0.0;
    }
    return regularized_gamma_p(df / 2.0, x / 2.0);
}

// Compute p-value from chi-squared statistic
auto chi_sq_pvalue(double stat, int df) -> double
{
    if (stat < 0.0 || std::isnan(stat) || std::isinf(stat))
    {
        return std::nan("");
    }
    return 1.0 - chi_squared_cdf(stat, df);
}

// Single effect Wald test (for additive-only or dominance-only)
auto wald_test_single(
    const Eigen::Ref<const Eigen::VectorXd>& z,
    const Eigen::Ref<const Eigen::VectorXd>& residual,
    const Eigen::Ref<const Eigen::MatrixXd>& v_inv) -> AssociationResult
{
    AssociationResult result;
    result.df = 1;

    // z'V⁻¹z
    const Eigen::VectorXd v_inv_z = v_inv * z;
    const double z_vinv_z = z.dot(v_inv_z);

    if (z_vinv_z < std::numeric_limits<double>::epsilon())
    {
        result.stat = 0.0;
        result.pvalue = 1.0;
        return result;
    }

    // β̂ = (z'V⁻¹z)⁻¹ z'V⁻¹r
    const double z_vinv_r = z.dot(v_inv * residual);
    const double beta = z_vinv_r / z_vinv_z;

    // SE(β̂) = sqrt((z'V⁻¹z)⁻¹)
    const double var_beta = 1.0 / z_vinv_z;
    const double se = std::sqrt(var_beta);

    // Wald statistic: β̂² / Var(β̂)
    result.stat = beta * beta / var_beta;
    result.pvalue = chi_sq_pvalue(result.stat, 1);

    return {
        .beta_a = beta,
        .se_a = se,
        .stat = result.stat,
        .pvalue = result.pvalue,
        .df = 1};
}

// Joint Wald test for additive + dominance (df = 2)
auto wald_test_joint(
    const Eigen::Ref<const Eigen::MatrixXd>& Z,
    const Eigen::Ref<const Eigen::VectorXd>& residual,
    const Eigen::Ref<const Eigen::MatrixXd>& v_inv) -> AssociationResult
{
    AssociationResult result;
    result.df = 2;

    // Z'V⁻¹Z (2×2)
    const Eigen::MatrixXd Z_vinv = v_inv * Z;
    const Eigen::Matrix2d Z_vinv_Z = Z.transpose() * Z_vinv;

    // Check for singularity
    const double det = Z_vinv_Z.determinant();
    if (std::abs(det) < std::numeric_limits<double>::epsilon())
    {
        result.stat = 0.0;
        result.pvalue = 1.0;
        return result;
    }

    // (Z'V⁻¹Z)⁻¹
    const Eigen::Matrix2d Z_vinv_Z_inv = Z_vinv_Z.inverse();

    // β̂ = (Z'V⁻¹Z)⁻¹ Z'V⁻¹r
    const Eigen::Vector2d Z_vinv_r = Z.transpose() * (v_inv * residual);
    const Eigen::Vector2d beta = Z_vinv_Z_inv * Z_vinv_r;

    // SE
    const double se_a = std::sqrt(Z_vinv_Z_inv(0, 0));
    const double se_d = std::sqrt(Z_vinv_Z_inv(1, 1));

    // Wald statistic: β̂' (Z'V⁻¹Z) β̂
    result.stat = beta.transpose() * Z_vinv_Z * beta;
    result.pvalue = chi_sq_pvalue(result.stat, 2);

    result.beta_a = beta(0);
    result.se_a = se_a;
    result.beta_d = beta(1);
    result.se_d = se_d;

    return result;
}

}  // namespace

auto wald_test(
    const EncodedSNP& snp,
    const Eigen::Ref<const Eigen::VectorXd>& residual,
    const Eigen::Ref<const Eigen::MatrixXd>& v_inv,
    AssocMode model) -> AssociationResult
{
    if (!snp.is_valid)
    {
        return AssociationResult{.stat = 0.0, .pvalue = 1.0, .df = 0};
    }

    switch (model)
    {
        case AssocMode::A:
        {
            auto result = wald_test_single(snp.Z.col(0), residual, v_inv);
            // beta_a is already set, keep beta_d as NaN
            return result;
        }
        case AssocMode::D:
        {
            auto result = wald_test_single(snp.Z.col(0), residual, v_inv);
            // Swap to beta_d since this is dominance model
            result.beta_d = result.beta_a;
            result.se_d = result.se_a;
            result.beta_a = std::nan("");
            result.se_a = std::nan("");
            return result;
        }
        case AssocMode::AD:
        {
            return wald_test_joint(snp.Z, residual, v_inv);
        }
    }

    return AssociationResult{};
}

}  // namespace gelex::gwas
