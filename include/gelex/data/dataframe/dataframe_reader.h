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
#include <format>
#include <fstream>
#include <numeric>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
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

struct ReadOptions
{
    std::vector<ColumnType> schema;
    char delimiter = '\t';
    bool header = true;
    std::vector<std::string> index_cols;
    std::vector<std::string> select_cols;
    std::unordered_set<
        std::string,
        TransparentHash<std::string>,
        TransparentEqual<std::string>>
        na_rep = {kDefaultNaRep.begin(), kDefaultNaRep.end()};
    NaAction na_action = NaAction::Exclude;
};

template <KeyType Key>
class DataFrameReader
{
   public:
    DataFrameReader(const std::string& path, const ReadOptions& options)
        : path_(path), options_(&options)
    {
    }
    auto read() -> DataFrame<Key>;

   private:
    auto parse_header() -> void;
    auto check_row() -> bool;
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
};

template <KeyType Key>
auto read_dataframe(const std::string& path, const ReadOptions& options)
    -> DataFrame<Key>;

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
        auto name = header_[select_pos_[i]];
        df.col_lookup_[name] = i;
        df.names_.push_back(name);
        df.columns_.emplace_back(name);
    }

    std::string line;
    while (std::getline(file, line))
    {
        tokenize(line);
        if (!check_row())
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
                    throw InvalidInputException(
                        "composite index requires string Key");
                }
                else
                {
                    df.index_.push_back(
                        std::format(
                            "{}{}{}",
                            tokens_[index_pos_[0]],
                            kSeparator,
                            tokens_[index_pos_[1]]));
                }
            }
        }

        for (std::size_t i = 0; i < select_pos_.size(); ++i)
        {
            parse_value(
                tokens_[select_pos_[i]], options_->schema[i], df.columns_[i]);
        }
    }

    if (options_->index_cols.empty())
    {
        if constexpr (std::is_same_v<Key, std::string>)
        {
            throw InvalidInputException(
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
        throw InvalidInputException(
            std::format("failed to parse '{}' as numeric", token));
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
auto DataFrameReader<Key>::parse_header() -> void
{
    std::unordered_map<std::string_view, std::size_t> name_to_pos;
    name_to_pos.reserve(header_.size());
    for (std::size_t i = 0; i < header_.size(); ++i)
    {
        name_to_pos[header_[i]] = i;
    }

    auto lookup = [&](const auto& col) -> std::size_t
    {
        auto it = name_to_pos.find(col);
        if (it == name_to_pos.end())
        {
            throw InvalidInputException(
                std::format("column not found: '{}'", col));
        }
        return it->second;
    };

    for (const auto& col : options_->index_cols)
    {
        index_pos_.push_back(lookup(col));
    }

    if (options_->select_cols.empty())
    {
        for (std::size_t i = 0; i < tokens_.size(); ++i)
        {
            if (std::ranges::find(index_pos_, i) == index_pos_.end())
            {
                select_pos_.push_back(i);
            }
        }
    }
    else
    {
        for (const auto& col : options_->select_cols)
        {
            select_pos_.push_back(lookup(col));
        }
    }

    assert(options_->schema.size() == select_pos_.size());
}

template <KeyType Key>
auto DataFrameReader<Key>::check_row() -> bool
{
    auto check = [this](std::size_t pos)
    {
        if (is_na(tokens_[pos]))
        {
            if (options_->na_action == NaAction::Throw)
            {
                throw InvalidInputException(
                    std::format("NA value found in column {}", pos));
            }
            return false;
        }
        return true;
    };
    for (auto pos : index_pos_)
    {
        if (!check(pos))
        {
            return false;
        }
    }
    for (auto pos : select_pos_)
    {
        if (!check(pos))
        {
            return false;
        }
    }
    return true;
}

template <KeyType Key>
auto read_dataframe(const std::string& path, const ReadOptions& options)
    -> DataFrame<Key>
{
    return DataFrameReader<Key>(path, options).read();
}

}  // namespace gelex::df

#endif  // GELEX_DATA_DATAFRAME_DATAFRAME_READER_H
