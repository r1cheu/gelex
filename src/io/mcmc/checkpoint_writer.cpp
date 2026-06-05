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

#include "gelex/io/mcmc/checkpoint_writer.h"

#include <array>
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
#include "gelex/io/binary_writer.h"

namespace gelex
{

namespace
{

class CheckpointWriter final : private infra::FieldVisitor
{
   public:
    explicit CheckpointWriter(io::BinaryWriter& writer) : writer_(writer) {}

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

    auto on(
        std::string_view name,
        Eigen::Ref<Eigen::VectorXi> value,
        FieldFlag flags) -> void override
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

    io::BinaryWriter& writer_;
};

}  // namespace

auto write_checkpoint(
    BayesState& state,
    const std::mt19937_64& rng,
    std::string_view prefix) -> void
{
    io::BinaryWriter writer(fmt::format("{}.ckpt", prefix));

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
