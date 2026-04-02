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

#ifndef GELEX_DATA_FRAME_DETAIL_DATAFRAME_READER_H_
#define GELEX_DATA_FRAME_DETAIL_DATAFRAME_READER_H_

#include <fmt/format.h>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/data/frame/dataframe_policy.h"
#include "gelex/data/frame/detail/row_validator.h"
#include "gelex/data/frame/detail/text_utils.h"
#include "gelex/data/frame/detail/value_parser.h"
#include "gelex/exception.h"
#include "gelex/types/sample_id.h"

namespace gelex
{

template <typename T>
class DataFrame;

namespace detail
{

template <typename T>
class DataFrameReader
{
   public:
    static constexpr size_t kFidColumnIndex = 0;
    static constexpr size_t kIidColumnIndex = 1;
    static constexpr size_t kFirstDataColumnIndex = 2;
    static constexpr size_t kRequiredInputColumns = 3;

    static auto load(
        const std::filesystem::path& path,
        const DataFrameLoadPolicy& policy) -> DataFrame<T>
    {
        DataFrameReader reader;
        return reader.load_impl(path, policy);
    }

   private:
    auto load_impl(
        const std::filesystem::path& path,
        const DataFrameLoadPolicy& policy) -> DataFrame<T>
    {
        std::ifstream file(path);
        if (!file.is_open())
        {
            throw GelexException(
                fmt::format("{}: failed to open file", path.string()));
        }

        std::string line;
        if (!std::getline(file, line))
        {
            throw GelexException(fmt::format("{}: is empty", path.string()));
        }
        if (!line.empty() && line.back() == '\r')
        {
            line.pop_back();
        }

        char delimiter = detect_delimiter(line);
        auto header_tokens = split_line_preserve_empty(line, delimiter);
        size_t expected_columns = header_tokens.size();
        size_t line_number = 1;

        validate_minimum_columns(
            expected_columns, kRequiredInputColumns, line_number);
        validate_header_prefix(header_tokens, line_number);

        auto all_names = extract_column_names(header_tokens);
        size_t num_data_columns = all_names.size();

        DataFrame<T> frame;
        if (policy.select_columns)
        {
            std::vector<std::string> selected_names;
            for (size_t col : *policy.select_columns)
            {
                if (col >= num_data_columns)
                {
                    throw GelexException(
                        fmt::format(
                            "select_columns index {} is out of range [0, {})",
                            col,
                            num_data_columns));
                }
                selected_names.push_back(std::move(all_names[col]));
                selected_input_columns_.push_back(col + kFirstDataColumnIndex);
            }
            frame.initialize_columns(std::move(selected_names));
        }
        else
        {
            frame.initialize_columns(std::move(all_names));
            for (size_t i = 0; i < num_data_columns; ++i)
            {
                selected_input_columns_.push_back(i + kFirstDataColumnIndex);
            }
        }
        row_buffer_.resize(frame.columns_.size());

        while (std::getline(file, line))
        {
            ++line_number;
            if (!line.empty() && line.back() == '\r')
            {
                line.pop_back();
            }
            auto tokens = split_line_preserve_empty(line, delimiter);
            append_row(frame, tokens, expected_columns, line_number, policy);
        }

        return frame;
    }

    auto append_row(
        DataFrame<T>& frame,
        std::span<const std::string_view> tokens,
        size_t expected_columns,
        size_t line_number,
        const DataFrameLoadPolicy& policy) -> void
    {
        validate_expected_columns(tokens.size(), expected_columns, line_number);
        validate_non_empty_token(tokens[kFidColumnIndex], line_number, "FID");
        validate_non_empty_token(tokens[kIidColumnIndex], line_number, "IID");

        if (!parse_data_columns(frame, tokens, line_number, policy))
        {
            return;
        }

        std::string fid
            = ValueParser<std::string>::parse(tokens[kFidColumnIndex]);
        std::string iid
            = ValueParser<std::string>::parse(tokens[kIidColumnIndex]);
        std::string index_key = make_sample_id(fid, iid);

        if (frame.index_map_.contains(index_key))
        {
            throw GelexException(
                fmt::format("line {} has duplicated index key", line_number));
        }

        for (size_t i = 0; i < frame.columns_.size(); ++i)
        {
            frame.columns_[i].append(row_buffer_[i]);
        }

        frame.index_map_.emplace(index_key, frame.index_.size());
        frame.index_.append(std::move(index_key));
    }

    static auto extract_column_names(
        std::span<const std::string_view> header_tokens)
        -> std::vector<std::string>
    {
        std::vector<std::string> names;
        names.reserve(header_tokens.size() - kFirstDataColumnIndex);
        for (size_t i = kFirstDataColumnIndex; i < header_tokens.size(); ++i)
        {
            names.emplace_back(header_tokens[i]);
        }
        return names;
    }

    static auto validate_header_prefix(
        std::span<const std::string_view> header_tokens,
        size_t line_number) -> void
    {
        if (header_tokens[kFidColumnIndex] != "FID"
            || header_tokens[kIidColumnIndex] != "IID")
        {
            throw GelexException(
                fmt::format(
                    "line {} header must start with FID and IID", line_number));
        }
    }

    auto parse_data_columns(
        const DataFrame<T>& frame,
        std::span<const std::string_view> tokens,
        size_t line_number,
        const DataFrameLoadPolicy& policy) -> bool
    {
        for (size_t i = 0; i < frame.columns_.size(); ++i)
        {
            size_t input_column = selected_input_columns_[i];
            size_t display_column = input_column - kFirstDataColumnIndex + 1;
            auto parsed_value = parse_data_token(
                tokens[input_column], line_number, display_column, policy);

            if (!parsed_value.has_value())
            {
                return false;
            }

            row_buffer_[i] = std::move(*parsed_value);
        }

        return true;
    }

    std::vector<T> row_buffer_;
    std::vector<size_t> selected_input_columns_;

    static auto parse_data_token(
        std::string_view token,
        size_t line_number,
        size_t column_index,
        const DataFrameLoadPolicy& policy) -> std::optional<T>
    {
        if (policy.missing_tokens.contains(token))
        {
            switch (policy.missing_value_action)
            {
                case MissingValueAction::Throw:
                    throw GelexException(
                        fmt::format(
                            "line {} column {} has missing value",
                            line_number,
                            column_index));
                case MissingValueAction::UseTypeDefault:
                    return T{};
                case MissingValueAction::SkipRow:
                    return std::nullopt;
            }
        }

        try
        {
            return ValueParser<T>::parse(token);
        }
        catch (const GelexException& ex)
        {
            throw GelexException(
                fmt::format(
                    "line {} column {} parse failed: {}",
                    line_number,
                    column_index,
                    ex.what()));
        }
    }
};

}  // namespace detail

}  // namespace gelex

#endif  // GELEX_DATA_FRAME_DETAIL_DATAFRAME_READER_H_
