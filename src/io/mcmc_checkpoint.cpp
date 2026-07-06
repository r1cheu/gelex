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

#include "gelex/io/mcmc_checkpoint.h"

#include <array>
#include <filesystem>
#include <random>
#include <span>
#include <sstream>
#include <string>
#include <string_view>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/bayes/state.h"
#include "gelex/exception.h"
#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/binary_writer.h"
#include "gelex/types/categorical_vector.h"

namespace gelex
{

namespace
{

class CheckpointReader final : private FieldVisitor
{
   public:
    explicit CheckpointReader(const BinaryReader& reader) : reader_(reader) {}

    auto read(BayesState& state) -> void { state.visit(*this); }

   private:
    auto on(
        std::string_view name,
        Eigen::Ref<Eigen::VectorXf> value,
        FieldFlag flags) -> void override
    {
        if (!has(flags, FieldFlag::checkpoint))
        {
            return;
        }

        const auto field = field_path(name);
        const auto stored = reader_.to_mat<float>(field);
        if (stored.rows() != value.rows() || stored.cols() != value.cols())
        {
            throw GelexException(
                fmt::format(
                    "CheckpointReader: field shape mismatch for {}: got {}x{}, "
                    "expected {}x{}",
                    field,
                    stored.rows(),
                    stored.cols(),
                    value.rows(),
                    value.cols()));
        }
        value = stored;
    }

    auto on(
        std::string_view name,
        Eigen::Ref<Eigen::VectorXd> value,
        FieldFlag flags) -> void override
    {
        if (!has(flags, FieldFlag::checkpoint))
        {
            return;
        }

        const auto field = field_path(name);
        const auto stored = reader_.to_mat<double>(field);
        if (stored.rows() != value.rows() || stored.cols() != value.cols())
        {
            throw GelexException(
                fmt::format(
                    "CheckpointReader: field shape mismatch for {}: got {}x{}, "
                    "expected {}x{}",
                    field,
                    stored.rows(),
                    stored.cols(),
                    value.rows(),
                    value.cols()));
        }
        value = stored;
    }

    auto on(std::string_view name, CategoricalVector& value, FieldFlag flags)
        -> void override
    {
        if (!has(flags, FieldFlag::checkpoint))
        {
            return;
        }

        const auto field = field_path(name);
        const auto stored = reader_.to_mat<int>(field);
        if (stored.rows() != value.rows() || stored.cols() != value.cols())
        {
            throw GelexException(
                fmt::format(
                    "CheckpointReader: field shape mismatch for {}: got {}x{}, "
                    "expected {}x{}",
                    field,
                    stored.rows(),
                    stored.cols(),
                    value.rows(),
                    value.cols()));
        }
        value = stored;
    }

    auto on(std::string_view name, double& value, FieldFlag flags)
        -> void override
    {
        if (!has(flags, FieldFlag::checkpoint))
        {
            return;
        }

        const auto field = field_path(name);
        const auto stored = reader_.to_mat<double>(field);
        if (stored.rows() != 1 || stored.cols() != 1)
        {
            throw GelexException(
                fmt::format(
                    "CheckpointReader: field shape mismatch for {}: got {}x{}, "
                    "expected 1x1",
                    field,
                    stored.rows(),
                    stored.cols()));
        }
        value = stored(0, 0);
    }

    auto on(std::string_view name, int& value, FieldFlag flags) -> void override
    {
        if (!has(flags, FieldFlag::checkpoint))
        {
            return;
        }

        const auto field = field_path(name);
        const auto stored = reader_.to_mat<int>(field);
        if (stored.rows() != 1 || stored.cols() != 1)
        {
            throw GelexException(
                fmt::format(
                    "CheckpointReader: field shape mismatch for {}: got {}x{}, "
                    "expected 1x1",
                    field,
                    stored.rows(),
                    stored.cols()));
        }
        value = stored(0, 0);
    }

    auto on(
        std::string_view name,
        std::span<const std::string>,
        FieldFlag flags) -> void override
    {
        if (!has(flags, FieldFlag::checkpoint))
        {
            return;
        }
        throw GelexException(
            "CheckpointReader: string list checkpoint field is not supported "
            "for "
            + field_path(name));
    }

    auto on(std::string_view name, std::string_view, FieldFlag flags)
        -> void override
    {
        if (!has(flags, FieldFlag::checkpoint))
        {
            return;
        }
        throw GelexException(
            "CheckpointReader: string checkpoint field is not supported for "
            + field_path(name));
    }

    const BinaryReader& reader_;
};

class CheckpointWriter final : private FieldVisitor
{
   public:
    explicit CheckpointWriter(BinaryWriter& writer) : writer_(writer) {}

    auto write(BayesState& state) -> void { state.visit(*this); }

   private:
    auto on(
        std::string_view name,
        Eigen::Ref<Eigen::VectorXf> value,
        FieldFlag flags) -> void override
    {
        if (!has(flags, FieldFlag::checkpoint))
        {
            return;
        }
        writer_.write(field_path(name), value);
    }

    auto on(
        std::string_view name,
        Eigen::Ref<Eigen::VectorXd> value,
        FieldFlag flags) -> void override
    {
        if (!has(flags, FieldFlag::checkpoint))
        {
            return;
        }
        writer_.write(field_path(name), value);
    }

    auto on(std::string_view name, CategoricalVector& value, FieldFlag flags)
        -> void override
    {
        if (!has(flags, FieldFlag::checkpoint))
        {
            return;
        }
        writer_.write(field_path(name), value);
    }

    auto on(std::string_view name, double& value, FieldFlag flags)
        -> void override
    {
        if (!has(flags, FieldFlag::checkpoint))
        {
            return;
        }
        writer_.write(field_path(name), value);
    }

    auto on(std::string_view name, int& value, FieldFlag flags) -> void override
    {
        if (!has(flags, FieldFlag::checkpoint))
        {
            return;
        }
        writer_.write(field_path(name), value);
    }

    auto on(
        std::string_view name,
        std::span<const std::string>,
        FieldFlag flags) -> void override
    {
        if (!has(flags, FieldFlag::checkpoint))
        {
            return;
        }
        throw GelexException(
            "CheckpointWriter: string list checkpoint field is not supported "
            "for "
            + field_path(name));
    }

    auto on(std::string_view name, std::string_view, FieldFlag flags)
        -> void override
    {
        if (!has(flags, FieldFlag::checkpoint))
        {
            return;
        }
        throw GelexException(
            "CheckpointWriter: string checkpoint field is not supported for "
            + field_path(name));
    }

    BinaryWriter& writer_;
};

}  // namespace

auto read_checkpoint(const std::filesystem::path& path, BayesState& state)
    -> std::mt19937_64
{
    BinaryReader reader(path.string());
    CheckpointReader checkpoint_reader{reader};
    checkpoint_reader.read(state);

    const auto rng_state = reader.to_strings("rng_state");
    if (rng_state.size() != 1)
    {
        throw GelexException(
            fmt::format(
                "{}: checkpoint rng_state count mismatch: got {}, expected 1",
                path.string(),
                rng_state.size()));
    }

    std::istringstream iss{std::string{rng_state.front()}};
    std::mt19937_64 rng{};
    iss >> rng;
    if (!iss)
    {
        throw GelexException(
            fmt::format(
                "{}: failed to decode checkpoint rng_state", path.string()));
    }
    return rng;
}

auto write_checkpoint(
    BayesState& state,
    const std::mt19937_64& rng,
    std::string_view prefix) -> void
{
    BinaryWriter writer(fmt::format("{}.ckpt", prefix));

    CheckpointWriter checkpoint_writer{writer};
    checkpoint_writer.write(state);

    std::ostringstream oss;
    oss << rng;
    const auto rng_state = oss.str();
    writer.write_strings(
        "rng_state",
        std::array<std::string_view, 1>{std::string_view{rng_state}});
    writer.close();
}

}  // namespace gelex
