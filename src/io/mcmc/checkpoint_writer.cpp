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

#include <fmt/format.h>
#include <cstdint>
#include <span>
#include <sstream>

#include "gelex/infra/record_visitor.h"
#include "gelex/io/detail/binary_writer.h"
#include "gelex/model/bayes/state.h"

namespace gelex
{

namespace
{

constexpr uint8_t kStateOnlyCheckpointVersion = 2;

auto write_scalar(
    io::detail::BinaryWriter& writer,
    std::string_view path,
    double value) -> void
{
    auto handle = writer.reserve<double>(path, 1, 1);
    writer.write(handle, value);
}

class CheckpointRecordSink final : public infra::RecordSink
{
   public:
    explicit CheckpointRecordSink(io::detail::BinaryWriter& writer)
        : writer_(writer)
    {
    }

    auto visit(
        std::string_view path,
        const Eigen::Ref<const Eigen::VectorXf>& value) -> void override
    {
        writer_.write(path, value);
    }

    auto visit(
        std::string_view path,
        const Eigen::Ref<const Eigen::VectorXd>& value) -> void override
    {
        writer_.write(path, value);
    }

    auto visit(
        std::string_view path,
        const Eigen::Ref<const Eigen::VectorXi>& value) -> void override
    {
        writer_.write(path, value);
    }

    auto visit(std::string_view path, const double& value) -> void override
    {
        write_scalar(writer_, path, value);
    }

   private:
    io::detail::BinaryWriter& writer_;
};

auto write_uint8(
    io::detail::BinaryWriter& writer,
    std::string_view path,
    uint8_t value) -> void
{
    auto handle = writer.reserve<uint8_t>(path, 1, 1);
    writer.write(handle, value);
}

auto write_rng(io::detail::BinaryWriter& writer, const std::mt19937_64& rng)
    -> void
{
    std::ostringstream oss;
    oss << rng;
    const auto str = oss.str();
    auto handle = writer.reserve<uint8_t>(
        "rng_state", static_cast<Eigen::Index>(str.size()), 1);
    writer.write(
        handle,
        std::span<const uint8_t>(
            reinterpret_cast<const uint8_t*>(str.data()), str.size()));
}

}  // namespace

auto write_checkpoint(
    const BayesState& state,
    const std::mt19937_64& rng,
    std::string_view prefix) -> void
{
    io::detail::BinaryWriter writer(fmt::format("{}.ckpt", prefix));

    write_uint8(writer, "format_version", kStateOnlyCheckpointVersion);
    CheckpointRecordSink sink(writer);
    state.visit_records(bayes::StateRecordSet::checkpoint, sink);
    write_rng(writer, rng);
}

}  // namespace gelex
