#pragma once

#include <algorithm>  // For std::ranges::sort, std::ranges::set_intersection
#include <filesystem>
#include <format>
#include <fstream>
#include <numeric>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>
//
#include <Eigen/Dense>

#include "../src/data/parser.h"

namespace gelex::detail
{

template <typename T>
concept CppStream = std::derived_from<T, std::ios_base>;

enum class file_type : uint8_t
{
    text,
    binary
};
template <CppStream StreamType>
[[nodiscard]] StreamType open_or_throw(
    const std::filesystem::path& path,
    file_type type = file_type::text)
{
    std::ios_base::openmode mode{};
    if constexpr (std::is_base_of_v<std::istream, StreamType>)
    {
        mode |= std::ios_base::in;
    }
    if constexpr (std::is_base_of_v<std::ostream, StreamType>)
    {
        mode |= std::ios_base::out;
        if constexpr (
            std::is_same_v<StreamType, std::ofstream>
            || std::is_same_v<StreamType, std::wofstream>)
        {
            mode |= std::ios_base::trunc;
        }
    }

    if (type == file_type::binary)
    {
        mode |= std::ios_base::binary;
    }

    StreamType stream(path, mode);

    if (!stream.is_open())
    {
        throw std::runtime_error(
            std::format("Failed to open file: '{}'", path.string()));
    }

    return stream;
}
void validate_path_or_throw(const std::filesystem::path& path);

std::unordered_set<std::string> get_ids_from_file(
    const std::filesystem::path& path,
    bool iid_only,
    bool skip_header = true);

}  // namespace gelex::detail
