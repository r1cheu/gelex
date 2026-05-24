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

#include "gelex/io/mcmc/checkpoint_reader.h"

#include <fmt/format.h>
#include <cstdint>
#include <random>
#include <sstream>
#include <string>

#include "gelex/infra/record_visitor.h"
#include "gelex/io/binary_reader.h"
#include "gelex/model/bayes/state.h"

namespace gelex
{

namespace
{

constexpr uint8_t kStateOnlyCheckpointVersion = 2;

auto read_scalar(const io::BinaryReader& reader, std::string_view path)
    -> double
{
    auto map = reader.to_map<double>(path);
    return map(0, 0);
}

auto read_uint8(const io::BinaryReader& reader, std::string_view path)
    -> uint8_t
{
    auto map = reader.to_map<uint8_t>(path);
    return map(0, 0);
}

class CheckpointRecordReader final : public infra::MutableRecordSink
{
   public:
    explicit CheckpointRecordReader(const io::BinaryReader& reader)
        : reader_(reader)
    {
    }

    auto visit(std::string_view path, Eigen::Ref<Eigen::VectorXf> value)
        -> void override
    {
        read_vector<float>(path, value);
    }

    auto visit(std::string_view path, Eigen::Ref<Eigen::VectorXd> value)
        -> void override
    {
        read_vector<double>(path, value);
    }

    auto visit(std::string_view path, Eigen::Ref<Eigen::VectorXi> value)
        -> void override
    {
        read_vector<int>(path, value);
    }

    auto visit(std::string_view path, double& value) -> void override
    {
        value = read_scalar(reader_, path);
    }

   private:
    template <typename T, typename Value>
    auto read_vector(std::string_view path, Eigen::Ref<Value> value) -> void
    {
        const auto stored = reader_.to_mat<T>(path);
        if (stored.cols() != 1 || stored.rows() != value.size())
        {
            throw GelexException(
                fmt::format(
                    "checkpoint record shape mismatch for {}: got {}x{}, "
                    "expected {}x1",
                    path,
                    stored.rows(),
                    stored.cols(),
                    value.size()));
        }
        value = stored.col(0);
    }

    const io::BinaryReader& reader_;
};

auto read_rng(const io::BinaryReader& reader) -> std::mt19937_64
{
    auto rng_map = reader.to_map<uint8_t>("rng_state");
    std::string str(
        reinterpret_cast<const char*>(rng_map.data()),
        static_cast<size_t>(rng_map.size()));
    std::istringstream iss(str);
    std::mt19937_64 rng{};
    iss >> rng;
    return rng;
}

}  // namespace

auto read_checkpoint(const std::filesystem::path& path, BayesState& state)
    -> std::mt19937_64
{
    io::BinaryReader reader(path.string());
    if (!reader.contains("format_version"))
    {
        throw GelexException(
            fmt::format(
                "{}: legacy checkpoint format is no longer supported for "
                "BayesState resume",
                path.string()));
    }

    const auto format_version = read_uint8(reader, "format_version");
    if (format_version != kStateOnlyCheckpointVersion)
    {
        throw GelexException(
            fmt::format(
                "{}: unsupported checkpoint format version {}",
                path.string(),
                format_version));
    }

    CheckpointRecordReader sink(reader);
    state.visit_records(bayes::StateRecordSet::checkpoint, sink);
    return read_rng(reader);
}

}  // namespace gelex
