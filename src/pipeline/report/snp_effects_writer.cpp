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

#include "gelex/pipeline/report/snp_effects_writer.h"

#include <format>
#include <memory>

#include "gelex/io/text_writer.h"

namespace gelex
{

using Eigen::Index;

SnpEffectsWriter::SnpEffectsWriter(
    const MCMCResult& result,
    const std::filesystem::path& bim_file_path,
    const std::filesystem::path& output_path)
    : result_(&result),
      bim_loader_(bim_file_path),
      writer_(std::make_unique<detail::TextWriter>(output_path))
{
}

SnpEffectsWriter::~SnpEffectsWriter() = default;

auto SnpEffectsWriter::write() -> void
{
    const auto* additive = result_->additive();
    if (additive == nullptr)
    {
        return;
    }

    write_header();

    for (Index i = 0; i < additive->coeffs.size(); ++i)
    {
        row_buf_.clear();
        write_snp_row(i);
        writer_->write(row_buf_);
    }
}

auto SnpEffectsWriter::write_header() -> void
{
    const auto* additive = result_->additive();
    const auto* dominant = result_->dominant();

    Index n_add_components = 0;
    Index n_dom_components = 0;
    if (additive != nullptr && additive->comp_probs.cols() > 0)
    {
        n_add_components = additive->comp_probs.cols();
    }
    if (dominant != nullptr && dominant->comp_probs.cols() > 0)
    {
        n_dom_components = dominant->comp_probs.cols();
    }

    std::string h
        = "Index\tID\tChrom\tPosition\tA1\tA2\tA1Freq"
          "\tAdd\tAddSE\tAddPVE";

    if (n_add_components > 2)
    {
        for (Index comp = 0; comp < n_add_components; ++comp)
        {
            h += std::format("\tpi_{}", comp);
        }
    }
    h += "\tPIP";

    if (dominant != nullptr)
    {
        h += "\tDom\tDomSE\tDomPVE";
        if (n_dom_components > 2)
        {
            for (Index comp = 0; comp < n_dom_components; ++comp)
            {
                h += std::format("\tpi_{}", comp);
            }
        }
        h += "\tPIP";
    }

    writer_->write(h);
}

auto SnpEffectsWriter::write_snp_row(Index snp_index) -> void
{
    const auto* additive = result_->additive();
    const auto* dominant = result_->dominant();

    row_buf_ += std::format("{}\t", snp_index + 1);

    write_snp_basic_info(snp_index);
    write_effects(additive, snp_index);
    write_component_probs(additive, snp_index);
    write_pip(additive, snp_index);
    write_effects(dominant, snp_index);
    write_component_probs(dominant, snp_index);
    write_pip(dominant, snp_index);
}

auto SnpEffectsWriter::write_snp_basic_info(Index snp_index) -> void
{
    if (snp_index < static_cast<Index>(bim_loader_.size()))
    {
        const auto& snp_info = bim_loader_.info()[snp_index];
        row_buf_ += std::format(
            "{}\t{}\t{}\t{}\t{}",
            snp_info.id,
            snp_info.chrom,
            snp_info.pos,
            snp_info.A1,
            snp_info.A2);

        if (result_->p_freq.size() > snp_index)
        {
            row_buf_ += std::format("\t{:.6f}", result_->p_freq(snp_index));
        }
        else
        {
            row_buf_ += "\tNA";
        }
    }
    else
    {
        row_buf_ += "\tNA\tNA\tNA\tNA\tNA";
    }
}

auto SnpEffectsWriter::write_effects(
    const BaseMarkerSummary* effect,
    Index snp_index) -> void
{
    if (effect == nullptr || snp_index >= effect->coeffs.size())
    {
        return;
    }

    row_buf_ += std::format(
        "\t{:.6f}\t{:.6f}",
        effect->coeffs.mean(snp_index),
        effect->coeffs.stddev(snp_index));

    if (effect->pve.size() > snp_index)
    {
        row_buf_ += std::format("\t{:.6e}", effect->pve.mean(snp_index));
    }
    else
    {
        row_buf_ += "\t0.0";
    }
}

auto SnpEffectsWriter::write_component_probs(
    const BaseMarkerSummary* effect,
    Index snp_index) -> void
{
    if (effect == nullptr || effect->comp_probs.cols() <= 2
        || snp_index >= effect->comp_probs.rows())
    {
        return;
    }

    for (Index comp = 0; comp < effect->comp_probs.cols(); ++comp)
    {
        row_buf_
            += std::format("\t{:.6f}", effect->comp_probs(snp_index, comp));
    }
}

auto SnpEffectsWriter::write_pip(
    const BaseMarkerSummary* effect,
    Index snp_index) -> void
{
    if (effect == nullptr)
    {
        return;
    }

    if (effect->pip.size() > snp_index)
    {
        row_buf_ += std::format("\t{:.6f}", effect->pip(snp_index));
    }
    else
    {
        row_buf_ += "\t1.0";
    }
}

}  // namespace gelex
