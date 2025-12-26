#include "predict/snp_predictor.h"

#include <cmath>
#include <format>

#include <omp.h>

#include "gelex/exception.h"

namespace gelex
{

SnpPredictor::SnpPredictor(const SnpEffects& effects) : effects_(effects) {}

SnpComputeResult SnpPredictor::compute_split(const Eigen::MatrixXd& genotype)
{
    const Eigen::Index n_samples = genotype.rows();
    const Eigen::Index n_snps = genotype.cols();

    Eigen::VectorXd frequencies = effects_.frequencies();
    Eigen::VectorXd additive_effects = effects_.additive_effects();
    const bool has_dominant = effects_.dominance_effects().size() > 0;
    Eigen::VectorXd dominant_effects
        = has_dominant ? effects_.dominance_effects() : Eigen::VectorXd();

    if (frequencies.size() != n_snps || additive_effects.size() != n_snps)
    {
        throw InvalidInputException(
            std::format(
                "Dimension mismatch: genotype has {} SNPs, but frequencies has "
                "{} and additive effects has {}",
                n_snps,
                frequencies.size(),
                additive_effects.size()));
    }
    if (has_dominant && dominant_effects.size() != n_snps)
    {
        throw InvalidInputException(
            std::format(
                "Dimension mismatch: genotype has {} SNPs, but dominant "
                "effects has {}",
                n_snps,
                dominant_effects.size()));
    }

    const double eps = 1e-10;

    Eigen::VectorXd add_predictions = Eigen::VectorXd::Zero(n_samples);
    Eigen::VectorXd dom_predictions = Eigen::VectorXd::Zero(n_samples);

#pragma omp parallel
    {
        Eigen::VectorXd thread_add = Eigen::VectorXd::Zero(n_samples);
        Eigen::VectorXd thread_dom = Eigen::VectorXd::Zero(n_samples);

#pragma omp for schedule(dynamic) nowait
        for (Eigen::Index j = 0; j < n_snps; ++j)
        {
            const double p = frequencies[j];
            const double q = 1.0 - p;

            const double scale_add = std::sqrt(std::max(2.0 * p * q, eps));
            const double mean_add = 2.0 * p;
            Eigen::VectorXd std_add
                = (genotype.col(j).array() - mean_add) / scale_add;
            thread_add += std_add * additive_effects[j];

            if (has_dominant)
            {
                const double mean_dom = 2.0 * p * p;
                const double scale_dom = std::max(2.0 * p * q, eps);

                Eigen::VectorXd dom_transformed = genotype.col(j).unaryExpr(
                    [p](double x) -> double
                    {
                        if (x == 1.0)
                        {
                            return 2.0 * p;
                        }
                        if (x == 2.0)
                        {
                            return (4.0 * p) - 2.0;
                        }
                        return 0.0;
                    });

                Eigen::VectorXd std_dom
                    = (dom_transformed.array() - mean_dom) / scale_dom;
                thread_dom += std_dom * dominant_effects[j];
            }
        }

#pragma omp critical
        {
            add_predictions += thread_add;
            dom_predictions += thread_dom;
        }
    }

    return SnpComputeResult{
        .add = std::move(add_predictions), .dom = std::move(dom_predictions)};
}

Eigen::VectorXd SnpPredictor::compute(const Eigen::MatrixXd& genotype)
{
    auto result = compute_split(genotype);
    return result.add + result.dom;
}

}  // namespace gelex
