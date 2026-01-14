#include "gelex/data/grm_code_policy.h"

#include <cmath>

#include <omp.h>
#include <Eigen/Core>

namespace gelex::grm
{

namespace detail
{
auto additive_mean_center(Eigen::Ref<Eigen::MatrixXd> genotype) -> void
{
#pragma omp parallel for default(none) shared(genotype)
    for (Eigen::Index i = 0; i < genotype.cols(); ++i)
    {
        genotype.col(i).array() -= genotype.col(i).mean();
    }
}
}  // namespace detail

auto Su::operator()(Eigen::Ref<Eigen::MatrixXd> genotype, bool use_additive)
    const -> void
{
    if (use_additive)
    {
        detail::additive_mean_center(genotype);
    }
    else
    {
#pragma omp parallel for default(none) shared(genotype)
        for (Eigen::Index i = 0; i < genotype.cols(); ++i)
        {
            auto col = genotype.col(i);
            double pA = col.mean() / 2.0;
            col = (col.array() == 2).select(0, col);
            col.array() -= 2 * pA * (1 - pA);
        }
    }
}

auto Zeng::operator()(Eigen::Ref<Eigen::MatrixXd> genotype, bool use_additive)
    const -> void
{
    if (use_additive)
    {
        detail::additive_mean_center(genotype);
    }
    else
    {
#pragma omp parallel for default(none) shared(genotype)
        for (Eigen::Index i = 0; i < genotype.cols(); ++i)
        {
            auto col = genotype.col(i);
            double pA = col.mean() / 2.0;
            double pa = 1 - pA;

            col = col.unaryExpr(
                [pA, pa](double val)
                {
                    if (val == 2)
                        return -2 * pa * pa;
                    if (val == 1)
                        return 2 * pA * pa;
                    if (val == 0)
                        return -2 * pA * pA;
                    return val;
                });
        }
    }
}

auto Yang::operator()(Eigen::Ref<Eigen::MatrixXd> genotype, bool use_additive)
    const -> void
{
    if (use_additive)
    {
#pragma omp parallel for default(none) shared(genotype)
        for (Eigen::Index i = 0; i < genotype.cols(); ++i)
        {
            auto col = genotype.col(i);
            double pA = col.mean() / 2.0;
            double pa = 1 - pA;
            double denom = std::sqrt(2 * pA * pa);

            if (denom < EPSILON)
            {
                col.setZero();
                continue;
            }

            col.array() -= 2 * pA;
            col.array() /= denom;
        }
    }
    else
    {
#pragma omp parallel for default(none) shared(genotype)
        for (Eigen::Index i = 0; i < genotype.cols(); ++i)
        {
            auto col = genotype.col(i);
            double pA = col.mean() / 2.0;
            double pa = 1 - pA;
            double denom = 2 * pA * pa;

            if (denom < EPSILON)
            {
                col.setZero();
                continue;
            }

            col = col.unaryExpr(
                [pA, pa, denom](double val)
                {
                    if (val == 2)
                        return -2 * pa * pa / denom;
                    if (val == 1)
                        return 2 * pA * pa / denom;
                    if (val == 0)
                        return -2 * pA * pA / denom;
                    return val;
                });
        }
    }
}

auto Vitezica::operator()(
    Eigen::Ref<Eigen::MatrixXd> genotype,
    bool use_additive) const -> void
{
    if (use_additive)
    {
#pragma omp parallel for default(none) shared(genotype)
        for (Eigen::Index i = 0; i < genotype.cols(); ++i)
        {
            auto col = genotype.col(i);
            auto nAA = static_cast<double>((col.array() == 2).count());
            auto nAa = static_cast<double>((col.array() == 1).count());
            col.array() -= (nAa + (2 * nAA));
        }
    }
    else
    {
#pragma omp parallel for default(none) shared(genotype)
        for (Eigen::Index i = 0; i < genotype.cols(); ++i)
        {
            auto col = genotype.col(i);

            auto nAA = static_cast<double>((col.array() == 2).count());
            auto nAa = static_cast<double>((col.array() == 1).count());
            auto naa = static_cast<double>((col.array() == 0).count());

            double denom = (nAA + naa - std::pow(nAA - naa, 2));
            if (denom < EPSILON)
            {
                col.setZero();
                continue;
            }

            col = col.unaryExpr(
                [nAA, nAa, naa, denom](double val)
                {
                    if (val == 2)
                        return -2 * naa * nAa / denom;
                    if (val == 1)
                        return 4 * nAA * naa / denom;
                    if (val == 0)
                        return -2 * nAA * nAa / denom;
                    return val;
                });
        }
    }
}

}  // namespace gelex::grm
