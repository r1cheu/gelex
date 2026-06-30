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

#include "gelex/io/reml.h"

#include <cstddef>
#include <iterator>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/data/sample_id.h"
#include "gelex/exception.h"
#include "gelex/freq/design.h"
#include "gelex/freq/model.h"
#include "gelex/io/detail/text_writer.h"
#include "gelex/types/fixed_designs.h"

namespace gelex::reml
{

auto write_summary(
    const FreqModel& model,
    const FreqState& state,
    double loglike,
    std::string_view prefix) -> void
{
    io::detail::TextWriter writer(fmt::format("{}.summary", prefix));
    writer.write_header(
        {"term", "type", "estimate", "se", "ratio", "ratio_se"});

    std::vector<std::string> fixed_names;
    fixed_names.reserve(static_cast<std::size_t>(model.fixed().X.cols()));
    for (std::size_t i = 0; i < model.fixed().names.size(); ++i)
    {
        if (model.fixed().levels[i])
        {
            for (const auto& level : *model.fixed().levels[i])
            {
                fixed_names.emplace_back(model.fixed().names[i] + "_" + level);
            }
        }
        else
        {
            fixed_names.emplace_back(model.fixed().names[i]);
        }
    }

    for (Eigen::Index i = 0; i < state.fixed().coeffs.size(); ++i)
    {
        const auto row = static_cast<std::size_t>(i);
        const auto term = row < fixed_names.size() ? fixed_names[row]
                                                   : fmt::format("X{}", i);
        writer.write(
            fmt::format(
                "{}\tfixed\t{:.8e}\t{:.8e}\t-\t-",
                term,
                state.fixed().coeffs(i),
                state.fixed().se(i)));
    }

    for (std::size_t i = 0; i < state.random().size(); ++i)
    {
        const auto& random = state.random()[i];
        writer.write(
            fmt::format(
                "{}\tvariance\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}",
                model.random()[i].name,
                random.variance,
                random.variance_se,
                random.variance_ratio,
                random.variance_ratio_se));
    }

    writer.write(
        fmt::format(
            "{}\tvariance\t{:.8e}\t{:.8e}\t-\t-",
            "Residual",
            state.residual().variance,
            state.residual().variance_se));

    writer.write(fmt::format("logL\tmodelfit\t{:.8e}\t-\t-\t-", loglike));
}

auto write_effects(
    const FreqModel& model,
    const FreqState& state,
    std::span<const std::string> sample_ids,
    std::string_view prefix) -> void
{
    if (sample_ids.size() != static_cast<std::size_t>(model.num_individuals()))
    {
        throw GelexException(
            fmt::format(
                "write_effects: sample IDs size {} != model individuals {}",
                sample_ids.size(),
                model.num_individuals()));
    }

    std::vector<std::string> fixed_names;
    fixed_names.reserve(static_cast<std::size_t>(model.fixed().X.cols()));
    for (std::size_t i = 0; i < model.fixed().names.size(); ++i)
    {
        if (model.fixed().levels[i])
        {
            for (const auto& level : *model.fixed().levels[i])
            {
                fixed_names.push_back(model.fixed().names[i] + "_" + level);
            }
        }
        else
        {
            fixed_names.push_back(model.fixed().names[i]);
        }
    }

    for (auto i = static_cast<Eigen::Index>(fixed_names.size());
         i < model.fixed().X.cols();
         ++i)
    {
        fixed_names.push_back(fmt::format("X{}", i));
    }

    std::vector<std::string> random_names;
    std::vector<const Eigen::VectorXd*> random_blups;
    random_names.reserve(state.random().size());
    random_blups.reserve(state.random().size());
    for (std::size_t i = 0; i < state.random().size(); ++i)
    {
        random_names.push_back(model.random()[i].name);
        random_blups.push_back(&state.random()[i].blup);
    }

    std::string header = "FID\tIID";
    for (Eigen::Index i = 0; i < model.fixed().X.cols(); ++i)
    {
        header += fmt::format("\t{}", fixed_names[static_cast<std::size_t>(i)]);
    }
    for (const auto& name : random_names)
    {
        header += fmt::format("\t{}", name);
    }
    header += "\tTOTAL";

    io::detail::TextWriter writer(fmt::format("{}.effects", prefix));
    writer.write(header);

    std::string line;
    line.reserve(128 + ((fixed_names.size() + random_blups.size()) * 16));
    for (Eigen::Index i = 0; i < model.num_individuals(); ++i)
    {
        const auto row = static_cast<std::size_t>(i);
        auto [fid, iid] = split_sample_id(sample_ids[row]);
        double total = 0.0;

        line.clear();
        fmt::format_to(std::back_inserter(line), "{}\t{}", fid, iid);
        for (Eigen::Index j = 0; j < model.fixed().X.cols(); ++j)
        {
            const double fixed_effect
                = model.fixed().X(i, j) * state.fixed().coeffs(j);
            total += fixed_effect;
            fmt::format_to(std::back_inserter(line), "\t{:.8e}", fixed_effect);
        }
        for (const auto* blup : random_blups)
        {
            total += (*blup)(i);
            fmt::format_to(std::back_inserter(line), "\t{:.8e}", (*blup)(i));
        }
        fmt::format_to(std::back_inserter(line), "\t{:.8e}", total);
        writer.write(line);
    }
}

auto write_loco_summary(
    const std::vector<LocoRemlResult>& results,
    std::string_view prefix) -> void
{
    io::detail::TextWriter writer(fmt::format("{}.loco.summary", prefix));
    writer.write_header(
        {"chr",
         "term",
         "type",
         "estimate",
         "se",
         "ratio",
         "ratio_se",
         "loglike",
         "converged"});

    for (const auto& result : results)
    {
        const int converged = result.converged ? 1 : 0;
        for (const auto& component : result.random)
        {
            writer.write(
                fmt::format(
                    "{}\t{}\tvariance\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}\t{:.8e}"
                    "\t{}",
                    result.chr_name,
                    component.name,
                    component.variance,
                    component.variance_se,
                    component.variance_ratio,
                    component.variance_ratio_se,
                    result.loglike,
                    converged));
        }
        writer.write(
            fmt::format(
                "{}\tResidual\tvariance\t{:.8e}\t{:.8e}\t-\t-\t{:.8e}\t{}",
                result.chr_name,
                result.residual_variance,
                result.residual_variance_se,
                result.loglike,
                converged));
    }
}

}  // namespace gelex::reml
