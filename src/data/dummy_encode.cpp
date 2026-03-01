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

#include "gelex/data/frame/dummy_encode.h"

#include <format>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/infra/logger.h"

namespace gelex
{

namespace
{

struct ColumnMeta
{
    std::unordered_map<std::string, Eigen::Index> level_to_id;
    std::vector<std::string> levels;
    Eigen::Index dummy_start_col = 0;
    bool should_emit = false;
};

class DummyEncoder
{
   public:
    explicit DummyEncoder(const DataFrame<std::string>& frame)
        : frame_(&frame),
          n_rows_(static_cast<Eigen::Index>(frame_->nrows())),
          n_data_cols_(frame_->ncols()),
          col_meta_(n_data_cols_)
    {
    }

    auto encode() -> DiscreteCovariate
    {
        collect_column_metadata();
        initialize_matrix();
        encode_matrix();
        return DiscreteCovariate{
            .names = std::move(names_),
            .levels = std::move(levels_),
            .reference_levels = std::move(reference_levels_),
            .X = std::move(X_)};
    }

   private:
    auto collect_column_metadata() -> void
    {
        for (size_t col_idx = 0; col_idx < n_data_cols_; ++col_idx)
        {
            const auto& column = frame_->column(col_idx);
            const auto& column_values = column.data();
            validate_column_size(column.name(), column_values.size());

            auto& meta = col_meta_[col_idx];
            meta.level_to_id.reserve(column_values.size());

            for (const auto& value : column_values)
            {
                register_level(meta, value);
            }

            if (meta.levels.size() < 2)
            {
                warn_monomorphic_column(column.name());
                continue;
            }

            register_emitted_column(meta, column.name());
        }
    }

    static auto warn_monomorphic_column(std::string_view column_name) -> void
    {
        const auto message = std::format(
            "column '{}' is monomorphic and will be skipped in dummy encoding",
            column_name);

        auto logger = logging::get();
        if (logger)
        {
            logger->warn(message);
            return;
        }

        std::clog << "[warn] " << message << '\n';
    }

    auto validate_column_size(std::string_view column_name, size_t size) const
        -> void
    {
        if (size == frame_->nrows())
        {
            return;
        }

        throw InvalidOperationException(
            std::format("column '{}' size mismatch", column_name));
    }

    static auto register_level(ColumnMeta& meta, const std::string& value)
        -> void
    {
        const auto inserted
            = meta.level_to_id
                  .emplace(value, static_cast<Eigen::Index>(meta.levels.size()))
                  .second;
        if (inserted)
        {
            meta.levels.push_back(value);
        }
    }

    auto register_emitted_column(ColumnMeta& meta, std::string_view column_name)
        -> void
    {
        meta.should_emit = true;
        meta.dummy_start_col = total_dummy_cols_;
        total_dummy_cols_ += static_cast<Eigen::Index>(meta.levels.size() - 1);

        names_.emplace_back(column_name);
        levels_.push_back(meta.levels);
        reference_levels_.push_back(meta.levels.front());
    }

    auto initialize_matrix() -> void
    {
        X_.resize(n_rows_, total_dummy_cols_);
        X_.setZero();
    }

    auto encode_matrix() -> void
    {
        for (size_t col_idx = 0; col_idx < n_data_cols_; ++col_idx)
        {
            const auto& meta = col_meta_[col_idx];
            if (!meta.should_emit)
            {
                continue;
            }

            const auto& column_values = frame_->column(col_idx).data();
            for (size_t row_idx = 0; row_idx < column_values.size(); ++row_idx)
            {
                const auto level_id
                    = meta.level_to_id.at(column_values[row_idx]);
                if (level_id == 0)
                {
                    continue;
                }

                const auto dummy_col = meta.dummy_start_col + (level_id - 1);
                X_(static_cast<Eigen::Index>(row_idx), dummy_col) = 1.0;
            }
        }
    }

    const DataFrame<std::string>* frame_;
    Eigen::Index n_rows_;
    size_t n_data_cols_;

    std::vector<ColumnMeta> col_meta_;

    std::vector<std::string> names_;
    std::vector<std::vector<std::string>> levels_;
    std::vector<std::string> reference_levels_;
    Eigen::MatrixXd X_;
    Eigen::Index total_dummy_cols_ = 0;
};

}  // namespace

auto DummyEncode(const DataFrame<std::string>& frame) -> DiscreteCovariate
{
    return DummyEncoder(frame).encode();
}

}  // namespace gelex
