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

#include "gelex/freq/gwas/single_tester.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cstddef>
#include <optional>
#include <unsupported/Eigen/SpecialFunctions>

#include "gelex/data/encode/encoder.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"
#include "gelex/data/genotype_method.h"
#include "gelex/freq/gwas/assoc_tester.h"
#include "gelex/freq/reml/operators.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/var.h"

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

auto SingleTester::run(
    const LocusEncoder& encoder,
    Eigen::Index start,
    const GwasOperators& reml) -> TestResults
{
    const Eigen::Index cols = Z_.cols();
    const EncodingSpec spec{encoding_spec_from_method(mode_, method_)};

#pragma omp parallel for schedule(static)
    for (Eigen::Index c = 0; c < cols; ++c)
    {
        const LocusStats stats{encoder.count(start + c)};
        encoder.expand(
            start + c, encoder.encoding(start + c, stats, spec), Z_.col(c));
        freqs_(c) = stats.has_nonmissing() ? stats.A1freq() : 0.0;
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
    };
}

auto SingleTester::wald_test(
    Eigen::Ref<Eigen::MatrixXd> Z,
    Eigen::Ref<Eigen::MatrixXd> W,
    const GwasOperators& reml,
    AssocOutput& output) -> void
{
    output.zt_Pr.noalias() = Z.transpose() * reml.Py;
    W.noalias() = reml.P * Z;
    output.zt_Pz = (Z.transpose() * W).diagonal();

    output.beta = (output.zt_Pr.array() / output.zt_Pz.array());
    output.pve
        = matvar(Z * output.beta.asDiagonal(), VarNormType::Population).array()
          / reml.Vp;
    output.se = (1.0 / output.zt_Pz.array()).sqrt();
    output.stats = (output.beta.array() / output.se.array()).square();
    output.p_value = (output.stats.array() * 0.5).sqrt().erfc();
}

}  // namespace gelex
