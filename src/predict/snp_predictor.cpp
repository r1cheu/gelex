#include "predict/snp_predictor.h"

#include <cmath>
#include <format>

#include <omp.h>

#include "gelex/exception.h"

namespace gelex
{

using Eigen::Index;
using Eigen::MatrixXd;
using Eigen::Ref;
using Eigen::VectorXd;

namespace
{
constexpr double kEpsilon = 1e-10;

std::pair<double, double> get_additive_params(double maf)
{
    double mean = 2.0 * maf;
    double variance = 2.0 * maf * (1.0 - maf);
    double scale = std::sqrt(std::max(variance, kEpsilon));
    return {mean, scale};
}

std::pair<double, double> get_dominance_params(double maf)
{
    double mean = 2.0 * maf * maf;
    double scale = 2.0 * maf * (1.0 - maf);
    scale = std::max(scale, kEpsilon);
    return {mean, scale};
}

double encode_dominance_value(double genotype_val, double maf)
{
    int g = static_cast<int>(std::round(genotype_val));
    if (g == 0)
    {
        return 0.0;
    }
    if (g == 1)
    {
        return 2.0 * maf;
    }
    if (g == 2)
    {
        return (4.0 * maf) - 2.0;
    }
    return 0.0;
}

template <bool HasDominance>
void run_prediction(
    const Ref<const MatrixXd>& genotype,
    const Ref<const VectorXd>& freqs,
    const Ref<const VectorXd>& add_effects,
    const Ref<const VectorXd>& dom_effects,
    Ref<VectorXd> add_dest,
    Ref<VectorXd> dom_dest)
{
    const Eigen::Index n_snps = genotype.cols();

    for (Eigen::Index j = 0; j < n_snps; ++j)
    {
        double p = freqs[j];

        auto [mu_add, sigma_add] = get_additive_params(p);
        double beta_add = add_effects[j];

        if (std::abs(beta_add) > kEpsilon)
        {
            add_dest.array()
                += (genotype.col(j).array() - mu_add) / sigma_add * beta_add;
        }

        if constexpr (HasDominance)
        {
            double beta_dom = dom_effects[j];
            if (std::abs(beta_dom) > kEpsilon)
            {
                auto [mu_dom, sigma_dom] = get_dominance_params(p);
                auto dom_encoder = [p](double val)
                { return encode_dominance_value(val, p); };
                dom_dest.array()
                    += (genotype.col(j).unaryExpr(dom_encoder).array() - mu_dom)
                       / sigma_dom * beta_dom;
            }
        }
    }
}
}  // namespace

SnpPredictor::SnpPredictor(const SnpEffects& effects) : effects_(effects) {}

void SnpPredictor::validate_dimensions(
    const Ref<const MatrixXd>& genotype) const
{
    Eigen::Index n_snps_data = genotype.cols();
    Eigen::Index n_snps_model = effects_.frequencies().size();

    if (n_snps_data != n_snps_model)
    {
        throw gelex::InvalidInputException(
            std::format(
                "Dimension mismatch: genotype matrix has {} columns (SNPs), "
                "but model expects {}.",
                n_snps_data,
                n_snps_model));
    }
}

SnpComputeResult SnpPredictor::compute(
    const Ref<const MatrixXd>& genotype) const
{
    validate_dimensions(genotype);

    const Eigen::Index n_samples = genotype.rows();

    bool has_dominance = (effects_.dominance_effects().size() > 0);

    Eigen::VectorXd add_scores = Eigen::VectorXd::Zero(n_samples);
    Eigen::VectorXd dom_scores;

    if (has_dominance)
    {
        dom_scores = Eigen::VectorXd::Zero(n_samples);
        run_prediction<true>(
            genotype,
            effects_.frequencies(),
            effects_.additive_effects(),
            effects_.dominance_effects(),
            add_scores,
            dom_scores);
        return SnpComputeResult{
            .add = std::move(add_scores), .dom = std::move(dom_scores)};
    }

    run_prediction<false>(
        genotype,
        effects_.frequencies(),
        effects_.additive_effects(),
        effects_.dominance_effects(),
        add_scores,
        dom_scores);

    return SnpComputeResult{
        .add = std::move(add_scores), .dom = std::move(dom_scores)};
}

}  // namespace gelex
