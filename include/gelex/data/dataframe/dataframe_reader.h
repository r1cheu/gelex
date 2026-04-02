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

#ifndef GELEX_DATA_DATAFRAME_DATAFRAME_READER_H
#define GELEX_DATA_DATAFRAME_DATAFRAME_READER_H

#include <algorithm>
#include <cassert>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <optional>
#include <ranges>
#include <string>
#include <string_view>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <type_traits>
#include <variant>
#include <vector>

#include "gelex/data/dataframe/column.h"
#include "gelex/data/dataframe/constants.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/key_type.h"
#include "gelex/data/dataframe/string_hash.h"
#include "gelex/io/parser.h"

namespace gelex::df
{
enum class NaAction : std::uint8_t
{
    Throw,
    Exclude
};

using Schema = std::variant<ColumnType, std::vector<ColumnType>>;

struct ReadOptions
{
    std::optional<Schema> schema;
    char delimiter = '\t';
    bool header = true;
    std::vector<std::size_t> index_cols;
    std::vector<std::size_t> select_cols;
    std::vector<std::string> names;
    StringSet na_rep = {kDefaultNaRep.begin(), kDefaultNaRep.end()};
    NaAction na_action = NaAction::Throw;
};

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
    auto read() -> DataFrame<Key>;

   private:
    auto parse_header() -> void;
    auto check_col_range(std::size_t idx) -> void;
    auto filter_row() -> bool;
    auto prepare_header(std::ifstream& file) -> void;

    template <typename T>
    static auto parse_arithmetic(std::string_view token) -> T;
    auto parse_value(std::string_view token, ColumnType type, Column& col)
        -> void;
    auto parse_key(std::string_view token) -> Key;
    auto is_na(std::string_view token) -> bool;
    auto tokenize(std::string_view line) -> void;

    std::filesystem::path path_;
    const ReadOptions* options_;
    std::size_t n_cols_ = 0;
    std::vector<std::string_view> tokens_;
    std::vector<std::string> header_;
    std::vector<std::size_t> select_pos_;
    std::vector<std::size_t> index_pos_;
    std::vector<ColumnType> resolved_schema_;
};

template <KeyType Key>
auto read_dataframe(
    const std::filesystem::path& path,
    const ReadOptions& options) -> DataFrame<Key>;

// --- Implementation ---

template <KeyType Key>
auto DataFrameReader<Key>::read() -> DataFrame<Key>
{
    auto file = detail::open_file<std::ifstream>(path_, std::ios::in);
    prepare_header(file);
    parse_header();

    DataFrame<Key> df;
    for (std::size_t i = 0; i < select_pos_.size(); ++i)
    {
        df.names_.emplace_back();
        df.columns_.emplace_back(header_[select_pos_[i]]);
        df.set_name(i, header_[select_pos_[i]]);
    }

    if (!options_->names.empty())
    {
        df.rename(options_->names);
    }

    std::string line;
    while (std::getline(file, line))
    {
        tokenize(line);
        if (!filter_row())
        {
            continue;
        }

        if (!options_->index_cols.empty())
        {
            if (index_pos_.size() == 1)
            {
                df.index_.push_back(parse_key(tokens_[index_pos_[0]]));
            }
            else
            {
                if constexpr (!std::is_same_v<Key, std::string>)
                {
                    throw GelexException("composite index requires string Key");
                }
                else
                {
                    df.index_.push_back(
                        fmt::format(
                            "{}",
                            fmt::join(
                                index_pos_
                                    | std::views::transform(
                                        [&](auto i) { return tokens_[i]; }),
                                std::string_view(&kSeparator, 1))));
                }
            }
        }

        for (std::size_t i = 0; i < select_pos_.size(); ++i)
        {
            parse_value(
                tokens_[select_pos_[i]], resolved_schema_[i], df.columns_[i]);
        }
    }

    if (options_->index_cols.empty())
    {
        if constexpr (std::is_same_v<Key, std::string>)
        {
            throw GelexException(
                "auto-generated index requires arithmetic Key");
        }
        else
        {
            auto n = df.columns_.empty() ? 0 : df.columns_[0].size();
            std::vector<Key> keys(n);
            std::iota(keys.begin(), keys.end(), Key{0});
            df.index_ = Index<Key>(std::move(keys));
        }
    }

    return df;
}

template <KeyType Key>
auto DataFrameReader<Key>::prepare_header(std::ifstream& file) -> void
{
    std::string line;
    std::getline(file, line);
    tokenize(line);
    n_cols_ = tokens_.size();
    if (!options_->header)
    {
        file.seekg(0);
        for (std::size_t i = 0; i < n_cols_; ++i)
        {
            header_.push_back(std::to_string(i));
        }
    }
    else
    {
        header_.reserve(tokens_.size());
        for (const auto& token : tokens_)
        {
            header_.emplace_back(token);
        }
    }
}

template <KeyType Key>
auto DataFrameReader<Key>::tokenize(std::string_view line) -> void
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
auto DataFrameReader<Key>::is_na(std::string_view token) -> bool
{
    return options_->na_rep.contains(token);
}

template <KeyType Key>
template <typename T>
auto DataFrameReader<Key>::parse_arithmetic(std::string_view token) -> T
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
auto DataFrameReader<Key>::parse_value(
    std::string_view token,
    ColumnType type,
    Column& col) -> void
{
    switch (type)
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
}

template <KeyType Key>
auto DataFrameReader<Key>::parse_key(std::string_view token) -> Key
{
    if constexpr (std::is_same_v<Key, std::string>)
    {
        return std::string(token);
    }
    else
    {
        return parse_arithmetic<Key>(token);
    }
}

template <KeyType Key>
auto DataFrameReader<Key>::check_col_range(std::size_t idx) -> void
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
auto DataFrameReader<Key>::parse_header() -> void
{
    for (auto idx : options_->index_cols)
    {
        check_col_range(idx);
        index_pos_.push_back(idx);
    }

    if (!options_->schema && options_->select_cols.empty())
    {
        return;
    }

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

    if (options_->schema)
    {
        std::visit(
            [this](const auto& s)
            {
                if constexpr (std::is_same_v<
                                  std::remove_cvref_t<decltype(s)>,
                                  ColumnType>)
                {
                    resolved_schema_.assign(select_pos_.size(), s);
                }
                else
                {
                    assert(
                        s.size() == select_pos_.size()
                        && "schema size must match selected columns");
                    resolved_schema_ = s;
                }
            },
            *options_->schema);
    }
}

// filter row based on NA values in index and selected columns
template <KeyType Key>
auto DataFrameReader<Key>::filter_row() -> bool
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
auto read_dataframe(
    const std::filesystem::path& path,
    const ReadOptions& options) -> DataFrame<Key>
{
    return DataFrameReader<Key>(path, options).read();
}

}  // namespace gelex::df

#endif  // GELEX_DATA_DATAFRAME_DATAFRAME_READER_H
