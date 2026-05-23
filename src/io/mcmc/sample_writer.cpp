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

#include "gelex/io/mcmc/sample_writer.h"

#include <fmt/format.h>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <variant>

#include <Eigen/Core>

#include "gelex/infra/record_visitor.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

namespace
{

using RecordHandle = mcmc::Writer::RecordHandle;

auto mode_name(GeneticMode mode) -> std::string_view
{
    for (const auto& [value, name] : kGeneticModeNames)
    {
        if (value == mode)
        {
            return name;
        }
    }
    return "unknown";
}

auto collect_mode_names(std::span<const GeneticMode> modes)
    -> std::vector<std::string_view>
{
    std::vector<std::string_view> names;
    names.reserve(modes.size());
    for (const auto mode : modes)
    {
        names.push_back(mode_name(mode));
    }
    return names;
}

auto write_sample_metadata(
    io::detail::BinaryWriter& writer,
    const BayesState& state) -> void
{
    std::vector<std::string_view> genetic_modes;
    genetic_modes.reserve(state.genetics().size());
    for (const auto& genetic : state.genetics())
    {
        genetic_modes.push_back(mode_name(genetic.type));
    }
    if (!genetic_modes.empty())
    {
        writer.write_strings("genetic/modes", genetic_modes);
    }

    for (const auto& [i, block] : std::views::enumerate(state.genetic_blocks()))
    {
        auto modes = collect_mode_names(block.modes());
        writer.write_strings(fmt::format("genetic_block/{}/modes", i), modes);
    }
}

class ReserveRecordSink final : public infra::RecordSink
{
   public:
    ReserveRecordSink(
        io::detail::BinaryWriter& writer,
        std::unordered_map<std::string, RecordHandle>& handles,
        Eigen::Index n_records)
        : writer_(writer), handles_(handles), n_records_(n_records)
    {
    }

    auto visit(
        std::string_view path,
        const Eigen::Ref<const Eigen::VectorXf>& value) -> void override
    {
        reserve<float>(path, value.size());
    }

    auto visit(
        std::string_view path,
        const Eigen::Ref<const Eigen::VectorXd>& value) -> void override
    {
        reserve<double>(path, value.size());
    }

    auto visit(
        std::string_view path,
        const Eigen::Ref<const Eigen::VectorXi>& value) -> void override
    {
        reserve<int>(path, value.size());
    }

    auto visit(std::string_view path, const double&) -> void override
    {
        reserve<double>(path, 1);
    }

   private:
    template <typename T>
    auto reserve(std::string_view path, Eigen::Index rows) -> void
    {
        handles_.emplace(
            std::string(path), writer_.reserve<T>(path, rows, n_records_));
    }

    io::detail::BinaryWriter& writer_;
    std::unordered_map<std::string, RecordHandle>& handles_;
    Eigen::Index n_records_;
};

class WriteRecordSink final : public infra::RecordSink
{
   public:
    WriteRecordSink(
        io::detail::BinaryWriter& writer,
        const std::unordered_map<std::string, RecordHandle>& handles)
        : writer_(writer), handles_(handles)
    {
    }

    auto visit(
        std::string_view path,
        const Eigen::Ref<const Eigen::VectorXf>& value) -> void override
    {
        write(path, value);
    }

    auto visit(
        std::string_view path,
        const Eigen::Ref<const Eigen::VectorXd>& value) -> void override
    {
        write(path, value);
    }

    auto visit(
        std::string_view path,
        const Eigen::Ref<const Eigen::VectorXi>& value) -> void override
    {
        write(path, value);
    }

    auto visit(std::string_view path, const double& value) -> void override
    {
        const auto& handle = require_handle<double>(path);
        writer_.write(handle, value);
    }

   private:
    template <typename Vector>
    auto write(std::string_view path, const Vector& value) -> void
    {
        using T = std::decay_t<decltype(value(0))>;
        const auto& handle = require_handle<T>(path);
        writer_.write(handle, value);
    }

    template <typename T>
    auto require_handle(std::string_view path) const
        -> const io::detail::SectionHandle<T>&
    {
        const auto it = handles_.find(std::string(path));
        if (it == handles_.end())
        {
            throw GelexException(
                fmt::format("sample record handle missing: {}", path));
        }
        const auto* handle
            = std::get_if<io::detail::SectionHandle<T>>(&it->second);
        if (handle == nullptr)
        {
            throw GelexException(
                fmt::format("sample record handle type mismatch: {}", path));
        }
        return *handle;
    }

    io::detail::BinaryWriter& writer_;
    const std::unordered_map<std::string, RecordHandle>& handles_;
};

}  // namespace

mcmc::Writer::Writer(
    const BayesState& state,
    std::string_view prefix,
    Eigen::Index n_records)
    : writer_(fmt::format("{}.samples", prefix))
{
    write_sample_metadata(writer_, state);
    ReserveRecordSink sink(writer_, record_handles_, n_records);
    state.visit_records(bayes::StateRecordSet::sample, sink);
}

void mcmc::Writer::write(const BayesState& state)
{
    WriteRecordSink sink(writer_, record_handles_);
    state.visit_records(bayes::StateRecordSet::sample, sink);
}

}  // namespace gelex
