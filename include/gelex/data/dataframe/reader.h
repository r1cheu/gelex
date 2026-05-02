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

#ifndef GELEX_DATA_DATAFRAME_READER_H
#define GELEX_DATA_DATAFRAME_READER_H

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include "gelex/data/dataframe/column.h"
#include "gelex/data/dataframe/constants.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/key_type.h"
#include "gelex/infra/string_hash.h"
#include "gelex/io/detail/parser.h"

namespace gelex::dataframe
{
enum class NaAction : std::uint8_t
{
    Throw,
    Exclude
};

struct ReadOptions
{
    char delimiter = '\t';
    bool header = true;
    std::vector<std::size_t> index_cols;
    std::vector<std::size_t> select_cols;
    std::vector<std::string> names;
    infra::StringSet na_rep = {kDefaultNaRep.begin(), kDefaultNaRep.end()};
    NaAction na_action = NaAction::Throw;
};

namespace detail
{

template <KeyType Key>
class DataFrameReader
{
   public:
    DataFrameReader(
        const std::filesystem::path& path,
        const ReadOptions& options)
        : path_(path), options_(&options)
    {
    }

    template <ValueType Val>
    auto read() -> DataFrame<Key>;
    auto read(std::span<const ColumnType> schema) -> DataFrame<Key>;
    auto read_index() -> Index<Key>;

   private:
    auto prepare(bool index_only = false) -> std::ifstream;

    template <typename ParseFn>
    auto scan_lines(std::ifstream& file, const ParseFn& parse)
        -> std::vector<Key>;

    auto build_index(std::vector<Key> keys, std::size_t n_rows) -> Index<Key>;
    auto build_header(std::ifstream& file) -> void;
    auto check_col_range(std::size_t idx) -> void;
    auto tokenize(std::string_view line) -> void;
    auto is_na(std::string_view token) -> bool;
    auto check_row() -> bool;
    auto append_key(std::vector<Key>& keys) -> void;

    template <typename T>
    static auto parse_arithmetic(std::string_view token) -> T;

    std::filesystem::path path_;
    const ReadOptions* options_;
    std::size_t n_cols_ = 0;
    std::vector<std::string_view> tokens_;
    std::vector<std::string> header_;
    std::vector<std::size_t> select_pos_;
    std::vector<std::size_t> index_pos_;
};

}  // namespace detail

// compile-time typed: all selected columns share the same Val type
template <KeyType Key, ValueType Val>
auto read_dataframe(
    const std::filesystem::path& path,
    const ReadOptions& options) -> DataFrame<Key>;

// per-column schema: runtime dispatch per column
template <KeyType Key>
auto read_dataframe(
    const std::filesystem::path& path,
    const ReadOptions& options,
    std::span<const ColumnType> schema) -> DataFrame<Key>;

// index-only: returns Index directly
template <KeyType Key>
auto read_index(const std::filesystem::path& path, const ReadOptions& options)
    -> Index<Key>;

// --- Implementation ---

template <KeyType Key>
auto detail::DataFrameReader<Key>::build_header(std::ifstream& file) -> void
{
    std::string line;
    std::getline(file, line);
    tokenize(line);
    n_cols_ = tokens_.size();
    if (!options_->header)
    {
        file.seekg(0);
        header_.reserve(n_cols_);
        for (std::size_t i = 0; i < n_cols_; ++i)
        {
            header_.push_back(std::to_string(i));
        }
    }
    else
    {
        header_.assign(tokens_.begin(), tokens_.end());
    }
}

template <KeyType Key>
auto detail::DataFrameReader<Key>::check_col_range(std::size_t idx) -> void
{
    if (idx >= n_cols_)
    {
        throw GelexException(
            fmt::format(
                "{} column index {} out of range (file has {} columns)",
                path_.string(),
                idx,
                n_cols_));
    }
}

template <KeyType Key>
auto detail::DataFrameReader<Key>::tokenize(std::string_view line) -> void
{
    tokens_.clear();
    if (!line.empty() && line.back() == '\r')
    {
        line.remove_suffix(1);
    }
    std::size_t pos = 0;
    while (pos <= line.size())
    {
        auto sep = line.find(options_->delimiter, pos);
        if (sep == std::string_view::npos)
        {
            tokens_.push_back(line.substr(pos));
            break;
        }
        tokens_.push_back(line.substr(pos, sep - pos));
        pos = sep + 1;
    }
}

template <KeyType Key>
auto detail::DataFrameReader<Key>::is_na(std::string_view token) -> bool
{
    return options_->na_rep.contains(token);
}

template <KeyType Key>
auto detail::DataFrameReader<Key>::check_row() -> bool
{
    auto check = [this](const auto& positions) -> bool
    {
        for (auto pos : positions)
        {
            if (is_na(tokens_[pos]))
            {
                if (options_->na_action == NaAction::Throw)
                {
                    throw GelexException(
                        fmt::format("NA value found in column {}", pos));
                }
                return false;
            }
        }
        return true;
    };
    return check(index_pos_) && check(select_pos_);
}

template <KeyType Key>
auto detail::DataFrameReader<Key>::append_key(std::vector<Key>& keys) -> void
{
    if (index_pos_.size() == 1)
    {
        auto token = tokens_[index_pos_[0]];
        if constexpr (std::is_same_v<Key, std::string>)
        {
            keys.push_back(std::string(token));
        }
        else
        {
            keys.push_back(parse_arithmetic<Key>(token));
        }
    }
    else
    {
        if constexpr (!std::is_same_v<Key, std::string>)
        {
            throw GelexException("composite index requires string Key");
        }
        else
        {
            keys.push_back(
                fmt::format(
                    "{}",
                    fmt::join(
                        index_pos_
                            | std::views::transform([&](auto i)
                                                    { return tokens_[i]; }),
                        std::string_view(&kSeparator, 1))));
        }
    }
}

template <KeyType Key>
template <typename T>
auto detail::DataFrameReader<Key>::parse_arithmetic(std::string_view token) -> T
{
    T val{};
    const auto* end = token.data() + token.size();
    auto [ptr, ec] = std::from_chars(token.data(), end, val);
    if (ec != std::errc{} || ptr != end)
    {
        throw GelexException(
            fmt::format("failed to parse '{}' as numeric", token));
    }
    return val;
}

template <KeyType Key>
auto detail::DataFrameReader<Key>::prepare(bool index_only) -> std::ifstream
{
    auto file
        = ::gelex::io::detail::open_file<std::ifstream>(path_, std::ios::in);
    build_header(file);
    // resolve index
    for (auto idx : options_->index_cols)
    {
        check_col_range(idx);
        index_pos_.push_back(idx);
    }

    // resolve columns
    if (!index_only)
    {
        if (options_->select_cols.empty())
        {
            for (std::size_t i = 0; i < n_cols_; ++i)
            {
                if (std::ranges::find(index_pos_, i) == index_pos_.end())
                {
                    select_pos_.push_back(i);
                }
            }
        }
        else
        {
            for (auto idx : options_->select_cols)
            {
                check_col_range(idx);
                select_pos_.push_back(idx);
            }
        }
    }
    return file;
}

template <KeyType Key>
template <typename ParseFn>
auto detail::DataFrameReader<Key>::scan_lines(
    std::ifstream& file,
    const ParseFn& parse) -> std::vector<Key>
{
    std::vector<Key> keys;
    std::size_t line_number = options_->header ? 1 : 0;
    std::string line;
    while (std::getline(file, line))
    {
        ++line_number;
        tokenize(line);
        if (!check_row())
        {
            continue;
        }
        if (!index_pos_.empty())
        {
            append_key(keys);
        }
        for (std::size_t i = 0; i < select_pos_.size(); ++i)
        {
            try
            {
                parse(i);
            }
            catch (const GelexException& e)
            {
                throw GelexException(
                    fmt::format(
                        "{}:{}:{}: {}",
                        path_.string(),
                        line_number,
                        select_pos_[i] + 1,
                        e.what()));
            }
        }
    }
    return keys;
}

template <KeyType Key>
auto detail::DataFrameReader<Key>::build_index(
    std::vector<Key> keys,
    std::size_t n_rows) -> Index<Key>
{
    if (!options_->index_cols.empty())
    {
        return Index<Key>(std::move(keys));
    }
    if constexpr (std::is_same_v<Key, std::string>)
    {
        throw GelexException("auto-generated index requires arithmetic Key");
    }
    else
    {
        return Index<Key>(
            std::views::iota(Key{0}, static_cast<Key>(n_rows))
            | std::ranges::to<std::vector>());
    }
}

template <KeyType Key>
template <ValueType Val>
auto detail::DataFrameReader<Key>::read() -> DataFrame<Key>
{
    auto file = prepare();
    std::vector<std::vector<Val>> col_data(select_pos_.size());

    auto keys = scan_lines(
        file,
        [&](std::size_t i)
        {
            if constexpr (std::is_same_v<Val, std::string>)
            {
                col_data[i].push_back(std::string(tokens_[select_pos_[i]]));
            }
            else
            {
                col_data[i].push_back(
                    parse_arithmetic<Val>(tokens_[select_pos_[i]]));
            }
        });

    DataFrame<Key> df;
    df.index_ = build_index(
        std::move(keys), col_data.empty() ? 0 : col_data[0].size());

    for (std::size_t i = 0; i < select_pos_.size(); ++i)
    {
        df.push_back(Column(header_[select_pos_[i]], std::move(col_data[i])));
    }
    if (!options_->names.empty())
    {
        df.rename(options_->names);
    }
    return df;
}

template <KeyType Key>
auto detail::DataFrameReader<Key>::read(std::span<const ColumnType> schema)
    -> DataFrame<Key>
{
    auto file = prepare();
    if (schema.size() != select_pos_.size())
    {
        throw GelexException(
            fmt::format(
                "schema size {} does not match selected columns {}",
                schema.size(),
                select_pos_.size()));
    }

    std::vector<Column> columns;
    columns.reserve(select_pos_.size());
    for (auto pos : select_pos_)
    {
        columns.emplace_back(header_[pos]);
    }

    auto keys = scan_lines(
        file,
        [&](std::size_t i)
        {
            auto token = tokens_[select_pos_[i]];
            auto& col = columns[i];
            switch (schema[i])
            {
                case ColumnType::Int:
                    col.push_back(parse_arithmetic<std::int32_t>(token));
                    break;
                case ColumnType::Float:
                    col.push_back(parse_arithmetic<float>(token));
                    break;
                case ColumnType::Double:
                    col.push_back(parse_arithmetic<double>(token));
                    break;
                case ColumnType::String:
                    col.push_back(std::string(token));
                    break;
            }
        });

    DataFrame<Key> df;
    df.index_
        = build_index(std::move(keys), columns.empty() ? 0 : columns[0].size());
    for (auto& col : columns)
    {
        df.push_back(std::move(col));
    }
    if (!options_->names.empty())
    {
        df.rename(options_->names);
    }
    return df;
}

template <KeyType Key>
auto detail::DataFrameReader<Key>::read_index() -> Index<Key>
{
    if (options_->index_cols.empty())
    {
        throw GelexException("read_index requires at least one index_col");
    }
    auto file = prepare(/*index_only=*/true);
    auto keys = scan_lines(file, [](std::size_t) {});
    return Index<Key>(std::move(keys));
}

// ================================================================
// read_dataframe / read_index entry points
// ================================================================

template <KeyType Key, ValueType Val>
auto read_dataframe(
    const std::filesystem::path& path,
    const ReadOptions& options) -> DataFrame<Key>
{
    return detail::DataFrameReader<Key>(path, options).template read<Val>();
}

template <KeyType Key>
auto read_dataframe(
    const std::filesystem::path& path,
    const ReadOptions& options,
    std::span<const ColumnType> schema) -> DataFrame<Key>
{
    return detail::DataFrameReader<Key>(path, options).read(schema);
}

template <KeyType Key>
auto read_index(const std::filesystem::path& path, const ReadOptions& options)
    -> Index<Key>
{
    return detail::DataFrameReader<Key>(path, options).read_index();
}

}  // namespace gelex::dataframe

#endif  // GELEX_DATA_DATAFRAME_READER_H
