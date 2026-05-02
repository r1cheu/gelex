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

#include "gelex/data/genotype/bed_pipe.h"

#include <fmt/format.h>
#include <limits>
#include <memory>
#include <span>

#include "gelex/data/dataframe/index.h"

#include <omp.h>

#include "data/bed_pipe/metadata.h"
#include "data/bed_pipe/mmap_reader.h"
#include "data/bed_pipe/sample_projection.h"
#include "data/bed_pipe/variant_decoder.h"
#include "gelex/exception.h"

namespace gelex::genotype
{

namespace
{

auto validate_chunk_range(
    Eigen::Index start_col,
    Eigen::Index end_col,
    Eigen::Index max_cols) -> void
{
    if (start_col < 0 || end_col > max_cols || start_col >= end_col)
    {
        throw gelex::GelexException(
            fmt::format(
                "invalid chunk range: [{}, {}). Total SNPs: {}",
                start_col,
                end_col,
                max_cols));
    }
}

auto validate_target_buffer_shape(
    const Eigen::Ref<Eigen::MatrixXd>& target_buf,
    Eigen::Index expected_rows,
    Eigen::Index expected_cols) -> void
{
    if (target_buf.rows() != expected_rows
        || target_buf.cols() != expected_cols)
    {
        throw std::runtime_error(
            "BedPipe::load_chunk: target_buf dimension mismatch");
    }
}

}  // namespace

BedPipe::BedPipe(
    const std::filesystem::path& bed_prefix,
    const dataframe::Index<std::string>& sample_index)
{
    auto metadata = gelex::genotype::detail::load_bed_metadata(bed_prefix);

    projection_ = std::make_unique<gelex::genotype::detail::SampleProjection>(
        metadata.raw_ids, sample_index);

    bed_reader_ = std::make_unique<gelex::genotype::detail::BedMmapReader>(
        metadata.bed_path, metadata.num_raw_snps, metadata.bytes_per_variant);

    decoder_ = std::make_unique<gelex::genotype::detail::BedVariantDecoder>(
        metadata.num_raw_samples,
        metadata.bytes_per_variant,
        projection_->mapping(),
        projection_->is_dense());

    num_raw_snps_ = metadata.num_raw_snps;
    num_output_samples_ = projection_->num_output_samples();
}

BedPipe::BedPipe(BedPipe&&) noexcept = default;

BedPipe& BedPipe::operator=(BedPipe&&) noexcept = default;

BedPipe::~BedPipe() = default;

auto BedPipe::load() const -> Eigen::MatrixXd
{
    return load_chunk(0, num_snps());
}

auto BedPipe::load_chunk(Eigen::Index start_col, Eigen::Index end_col) const
    -> Eigen::MatrixXd
{
    const Eigen::Index max_cols = num_snps();
    validate_chunk_range(start_col, end_col, max_cols);

    const Eigen::Index rows = num_samples();
    const Eigen::Index cols = end_col - start_col;

    Eigen::MatrixXd result(rows, cols);
    load_chunk(result, start_col, end_col);
    return result;
}

void BedPipe::load_chunk(
    Eigen::Ref<Eigen::MatrixXd> target_buf,
    Eigen::Index start_col,
    Eigen::Index end_col) const
{
    const Eigen::Index max_cols = num_snps();
    validate_chunk_range(start_col, end_col, max_cols);

    const Eigen::Index num_output_rows = num_samples();
    const Eigen::Index num_output_cols = end_col - start_col;
    validate_target_buffer_shape(target_buf, num_output_rows, num_output_cols);

    if (!projection_->is_dense())
    {
        target_buf.setConstant(std::numeric_limits<double>::quiet_NaN());
    }

    const uint8_t* chunk_ptr
        = bed_reader_->chunk_ptr(start_col, num_output_cols);
    if (chunk_ptr == nullptr)
    {
        throw gelex::GelexException(
            "BedPipe::load_chunk: mapped BED payload is truncated");
    }

    const Eigen::Index bytes_per_variant = bed_reader_->bytes_per_variant();

#pragma omp parallel for schedule(static)
    for (Eigen::Index j = 0; j < num_output_cols; ++j)
    {
        const uint8_t* src_ptr
            = chunk_ptr + static_cast<size_t>(j * bytes_per_variant);

        double* col_data_ptr = target_buf.col(j).data();
        std::span<double> target_span(col_data_ptr, num_output_rows);
        decoder_->decode(src_ptr, target_span);
    }
}

auto BedPipe::num_samples() const -> Eigen::Index
{
    return num_output_samples_;
}

auto BedPipe::num_snps() const -> Eigen::Index
{
    return num_raw_snps_;
}

auto BedPipe::select(std::span<const Eigen::Index> col_indices) const
    -> Eigen::MatrixXd
{
    const Eigen::Index num_output_rows = num_samples();
    const auto num_output_cols = static_cast<Eigen::Index>(col_indices.size());

    Eigen::MatrixXd result(num_output_rows, num_output_cols);
    result.setConstant(std::numeric_limits<double>::quiet_NaN());

#pragma omp parallel for schedule(static)
    for (Eigen::Index j = 0; j < num_output_cols; ++j)
    {
        const Eigen::Index snp_idx = col_indices[static_cast<size_t>(j)];
        if (snp_idx < 0 || snp_idx >= num_raw_snps_)
        {
            continue;
        }

        const uint8_t* src_ptr = bed_reader_->chunk_ptr(snp_idx, 1);
        if (src_ptr == nullptr)
        {
            continue;
        }

        double* col_data_ptr = result.col(j).data();
        std::span<double> target_span(col_data_ptr, num_output_rows);
        decoder_->decode(src_ptr, target_span);
    }

    return result;
}

}  // namespace gelex::genotype
