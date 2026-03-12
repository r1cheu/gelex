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

#include "gelex/model/bayes/writer/mcmc_writer.h"

#include <cstdint>
#include <format>
#include <ranges>
#include <string_view>

#include <Eigen/Core>

#include "gelex/io/binary_format.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/states.h"

namespace gelex
{

using detail::binary_format::dtype_code;

MCMCWriter::MCMCWriter(
    const BayesModel& model,
    std::string_view prefix,
    Eigen::Index n_records)
    : writer_(std::format("{}.samples", prefix))
{
    const auto cols = static_cast<uint64_t>(n_records);

    // Fixed
    fixed_coeffs_ = writer_.reserve(
        {detail::EffectType::Fixed, detail::DataKind::Coeff},
        dtype_code<double>(),
        static_cast<uint64_t>(model.fixed().X.cols()),
        cols);

    // Random
    for (uint8_t i = 0; i < static_cast<uint8_t>(model.random().size()); ++i)
    {
        const auto n_coeffs = static_cast<uint64_t>(model.random()[i].X.cols());
        auto coeffs_h = writer_.reserve(
            {detail::EffectType::Random, detail::DataKind::Coeff, i},
            dtype_code<double>(),
            n_coeffs,
            cols);
        auto variance_h = writer_.reserve(
            {detail::EffectType::Random, detail::DataKind::Variance, i},
            dtype_code<double>(),
            1,
            cols);
        random_.push_back({.coeffs = coeffs_h, .variance = variance_h});
    }

    // Genetic
    for (const auto& effect : model.genetics())
    {
        auto sect = detail::to_section_effect_type(effect.type);
        const auto n_snps = static_cast<uint64_t>(bayes::get_cols(effect.X));

        GeneticHandles gh;
        gh.section_effect = sect;
        gh.coeffs = writer_.reserve(
            {sect, detail::DataKind::Coeff},
            dtype_code<double>(),
            n_snps,
            cols);
        gh.variance = writer_.reserve(
            {sect, detail::DataKind::Variance}, dtype_code<double>(), 1, cols);

        if (effect.mixture)
        {
            gh.mixture_tracker = writer_.reserve(
                {sect, detail::DataKind::Mixture},
                dtype_code<int8_t>(),
                n_snps,
                cols);
            if (effect.mixture->estimate_pi)
            {
                const auto n_pi = static_cast<uint64_t>(
                    effect.mixture->init_proportion.size());
                gh.pi = writer_.reserve(
                    {sect, detail::DataKind::Pi},
                    dtype_code<double>(),
                    n_pi,
                    cols);
            }
        }
        if (effect.sign)
        {
            gh.sign_tracker = writer_.reserve(
                {sect, detail::DataKind::Sign},
                dtype_code<int8_t>(),
                n_snps,
                cols);
            gh.positive_prob = writer_.reserve(
                {sect, detail::DataKind::PositiveProb},
                dtype_code<double>(),
                1,
                cols);
        }

        genetic_.push_back(gh);
    }

    // Residual
    residual_variance_ = writer_.reserve(
        {detail::EffectType::Residual, detail::DataKind::Variance},
        dtype_code<double>(),
        1,
        cols);
}

void MCMCWriter::write(const BayesState& state)
{
    // Fixed
    writer_.write(
        fixed_coeffs_,
        reinterpret_cast<const char*>(state.fixed().coeffs.data()),
        static_cast<std::streamsize>(
            state.fixed().coeffs.size() * sizeof(double)));

    // Random
    for (auto&& [handles, rs] : std::views::zip(random_, state.random()))
    {
        writer_.write(
            handles.coeffs,
            reinterpret_cast<const char*>(rs.coeffs.data()),
            static_cast<std::streamsize>(rs.coeffs.size() * sizeof(double)));
        writer_.write(
            handles.variance,
            reinterpret_cast<const char*>(&rs.variance),
            sizeof(double));
    }

    // Genetic
    for (auto&& [gh, gs] : std::views::zip(genetic_, state.genetics()))
    {
        writer_.write(
            gh.coeffs,
            reinterpret_cast<const char*>(gs.coeffs.data()),
            static_cast<std::streamsize>(gs.coeffs.size() * sizeof(double)));
        writer_.write(
            gh.variance,
            reinterpret_cast<const char*>(&gs.variance),
            sizeof(double));

        if (gh.mixture_tracker && gs.mixture)
        {
            writer_.write(
                *gh.mixture_tracker,
                reinterpret_cast<const char*>(gs.mixture->tracker.data()),
                static_cast<std::streamsize>(
                    gs.mixture->tracker.size() * sizeof(int8_t)));
            if (gh.pi)
            {
                writer_.write(
                    *gh.pi,
                    reinterpret_cast<const char*>(
                        gs.mixture->pi.proportion.data()),
                    static_cast<std::streamsize>(
                        gs.mixture->pi.proportion.size() * sizeof(double)));
            }
        }
        if (gh.sign_tracker && gs.sign)
        {
            writer_.write(
                *gh.sign_tracker,
                reinterpret_cast<const char*>(gs.sign->tracker.data()),
                static_cast<std::streamsize>(
                    gs.sign->tracker.size() * sizeof(int8_t)));
            writer_.write(
                *gh.positive_prob,
                reinterpret_cast<const char*>(&gs.sign->positive_prob),
                sizeof(double));
        }
    }

    // Residual
    writer_.write(
        residual_variance_,
        reinterpret_cast<const char*>(&state.residual().variance),
        sizeof(double));
}

void MCMCWriter::finalize()
{
    writer_.finalize();
}

}  // namespace gelex
