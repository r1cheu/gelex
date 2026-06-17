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

#include "gelex/algo/gwas/single_tester.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cstddef>
#include <optional>
#include <unsupported/Eigen/SpecialFunctions>

#include "gelex/algo/gwas/assoc_tester.h"
#include "gelex/algo/reml/result.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/locus_encoding.h"
#include "gelex/infra/stats/detail/var.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

SingleTester::SingleTester(GeneticMode mode, GenotypeMethod method)
    : mode_(mode), method_(method)
{
}

auto SingleTester::resize(Eigen::Index n_samples, Eigen::Index chunk_size)
    -> void
{
    Z_.resize(n_samples, chunk_size);
    W_.resize(n_samples, chunk_size);
    freqs_.resize(chunk_size);
    output_.resize(chunk_size);
}

auto SingleTester::genotype_buffer() -> Eigen::Ref<Eigen::MatrixXd>
{
    return Z_;
}

auto SingleTester::run(const RemlResult& reml) -> TestResults
{
    if (mode_ == GeneticMode::A)
    {
        const LociEncoding encoding{
            encode_inplace<double>(Z_, GeneticMode::A, method_)};
        for (const LocusEncoding& locus : encoding.loci)
        {
            freqs_(locus.column_index)
                = locus.stats.has_nonmissing() ? locus.stats.A1freq() : 0.0;
        }
    }
    else
    {
        const LociEncoding encoding{
            encode_inplace<double>(Z_, GeneticMode::D, method_)};
        for (const LocusEncoding& locus : encoding.loci)
        {
            freqs_(locus.column_index)
                = locus.stats.has_nonmissing() ? locus.stats.A1freq() : 0.0;
        }
    }

    wald_test(Z_, W_, reml, output_);

    const auto n = static_cast<size_t>(Z_.cols());
    return {
        .freq = {freqs_.data(), n},
        .additive = {
            .beta = {output_.beta.data(), n},
            .se = {output_.se.data(), n},
            .p = {output_.p_value.data(), n},
            .pve = {output_.pve.data(), n},
        },
        .dominance = std::nullopt,
        .joint_p = std::nullopt,
        .total_pve = std::nullopt,
        .joint_covariance = std::nullopt,
    };
}

auto SingleTester::wald_test(
    Eigen::Ref<Eigen::MatrixXd> Z,
    Eigen::Ref<Eigen::MatrixXd> W,
    const RemlResult& reml,
    AssocOutput& output) -> void
{
    output.zt_Pr.noalias() = Z.transpose() * reml.Py;
    W.noalias() = reml.P * Z;
    output.zt_Pz = (Z.transpose() * W).diagonal();

    output.beta = (output.zt_Pr.array() / output.zt_Pz.array());
    output.pve = stats::detail::matvar(
                     Z * output.beta.asDiagonal(),
                     stats::detail::VarNormType::Population)
                     .array()
                 / reml.Vp;
    output.se = (1.0 / output.zt_Pz.array()).sqrt();
    output.stats = (output.beta.array() / output.se.array()).square();
    output.p_value = (output.stats.array() * 0.5).sqrt().erfc();
}

}  // namespace gelex
