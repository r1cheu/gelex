/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "gelex/algo/gwas/joint_tester.h"

#include <cmath>
#include <limits>

#include <Eigen/Dense>
#include <unsupported/Eigen/SpecialFunctions>

#include "gelex/data/genotype/genotype_processor.h"
#include "gelex/infra/stats/descriptive.h"

namespace gelex
{

JointTester::JointTester(GenotypeProcessMethod method) : method_(method) {}

auto JointTester::resize(Eigen::Index n_samples, Eigen::Index chunk_size)
    -> void
{
    raw_.resize(n_samples, chunk_size);
    Z_a_.resize(n_samples, chunk_size);
    Z_d_.resize(n_samples, chunk_size);
    W_.resize(n_samples, chunk_size);
    freqs_.resize(chunk_size);

    add_.resize(chunk_size);
    dom_.resize(chunk_size);
    zt_a_Pzd_.resize(chunk_size);
    total_pve_.resize(chunk_size);
}

auto JointTester::genotype_buffer() -> Eigen::Ref<Eigen::MatrixXd>
{
    return raw_;
}

auto JointTester::run(const RemlResult& reml) -> TestResults
{
    Z_a_ = raw_;
    Z_d_ = raw_;
    gelex::genotype::process_matrix<GeneticMode::A>(method_, Z_a_, &freqs_);
    gelex::genotype::process_matrix<GeneticMode::D>(method_, Z_d_);

    // W = P * Z_a  →  cross-products involving Z_a
    W_.noalias() = reml.P * Z_a_;
    add_.zt_Pr.noalias() = Z_a_.transpose() * reml.Py;
    add_.zt_Pz = (Z_a_.transpose() * W_).diagonal();
    zt_a_Pzd_ = (Z_d_.transpose() * W_).diagonal();

    // W = P * Z_d  →  cross-products involving Z_d
    W_.noalias() = reml.P * Z_d_;
    dom_.zt_Pr.noalias() = Z_d_.transpose() * reml.Py;
    dom_.zt_Pz = (Z_d_.transpose() * W_).diagonal();

    // per-SNP 2×2 solve and Wald tests
    const auto n = Z_a_.cols();
    constexpr auto nan = std::numeric_limits<double>::quiet_NaN();

    for (Eigen::Index i = 0; i < n; ++i)
    {
        Eigen::Matrix2d XtPX;
        XtPX << add_.zt_Pz(i), zt_a_Pzd_(i), zt_a_Pzd_(i), dom_.zt_Pz(i);
        Eigen::Vector2d XtPy(add_.zt_Pr(i), dom_.zt_Pr(i));

        if (std::abs(XtPX.determinant()) < 1e-30)
        {
            add_.beta(i) = add_.se(i) = add_.p_value(i) = nan;
            dom_.beta(i) = dom_.se(i) = dom_.p_value(i) = nan;
            continue;
        }

        Eigen::Matrix2d inv = XtPX.inverse();
        Eigen::Vector2d beta = inv * XtPy;

        add_.beta(i) = beta(0);
        dom_.beta(i) = beta(1);

        add_.se(i) = std::sqrt(inv(0, 0));
        dom_.se(i) = std::sqrt(inv(1, 1));

        // df=1 Wald: stat = beta^2 / var(beta)
        const double stat_a = (beta(0) * beta(0)) / inv(0, 0);
        const double stat_d = (beta(1) * beta(1)) / inv(1, 1);
        add_.p_value(i) = std::erfc(std::sqrt(stat_a * 0.5));
        dom_.p_value(i) = std::erfc(std::sqrt(stat_d * 0.5));
    }

    // per-effect PVE
    add_.pve
        = stats::detail::var(Z_a_ * add_.beta.asDiagonal()).transpose().array()
          / reml.Vp;
    dom_.pve
        = stats::detail::var(Z_d_ * dom_.beta.asDiagonal()).transpose().array()
          / reml.Vp;

    // total PVE: var(Z_a * beta_a + Z_d * beta_d) / Vp
    Eigen::MatrixXd pred
        = Z_a_ * add_.beta.asDiagonal() + Z_d_ * dom_.beta.asDiagonal();
    total_pve_ = stats::detail::var(pred).transpose().array() / reml.Vp;

    const auto n_snps = static_cast<size_t>(n);
    return {
        .freq = {freqs_.data(), n_snps},
        .additive = {
            .beta = {add_.beta.data(), n_snps},
            .se = {add_.se.data(), n_snps},
            .p = {add_.p_value.data(), n_snps},
            .pve = {add_.pve.data(), n_snps},
        },
        .dominance = TestResult{
            .beta = {dom_.beta.data(), n_snps},
            .se = {dom_.se.data(), n_snps},
            .p = {dom_.p_value.data(), n_snps},
            .pve = {dom_.pve.data(), n_snps},
        },
        .total_pve = std::span{total_pve_.data(), n_snps},
    };
}

}  // namespace gelex
