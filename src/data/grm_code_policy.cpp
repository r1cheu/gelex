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

            col = (col.array() == 2).select(-2 * pa * pa, col);
            col = (col.array() == 1).select(2 * pA * pa, col);
            col = (col.array() == 0).select(-2 * pA * pA, col);
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

            col = (col.array() == 2).select(-2 * pa * pa, col);
            col = (col.array() == 1).select(2 * pA * pa, col);
            col = (col.array() == 0).select(-2 * pA * pA, col);
            col.array() /= denom;
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

            double nAA = 0;
            double nAa = 0;
            double naa = 0;
            for (Eigen::Index j = 0; j < col.size(); ++j)
            {
                const double val = col(j);
                if (val == 2)
                {
                    ++nAA;
                }
                else if (val == 1)
                {
                    ++nAa;
                }
                else if (val == 0)
                {
                    ++naa;
                }
            }

            double denom = (nAA + naa - std::pow(nAA - naa, 2));
            if (denom < EPSILON)
            {
                col.setZero();
                continue;
            }

            col = (col.array() == 2).select(-2 * naa * nAa, col);
            col = (col.array() == 1).select(4 * nAA * naa, col);
            col = (col.array() == 0).select(-2 * nAA * nAa, col);
            col.array() /= denom;
        }
    }
}

}  // namespace gelex::grm
