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

#include "gelex/model/bayes/writer/snp_effects_writer.h"

#include <cstdint>
#include <memory>
#include <string>

#include <fmt/format.h>

#include "gelex/data/reader.h"
#include "gelex/io/detail/text_writer.h"

namespace gelex
{

using Eigen::Index;

SnpEffectsWriter::SnpEffectsWriter(
    const mcmc::Result& result,
    const std::filesystem::path& bim_file_path,
    const std::filesystem::path& output_path)
    : result_(&result),
      bim_(read_bim(bim_file_path)),
      writer_(std::make_unique<io::detail::TextWriter>(output_path))
{
    bim_keys_ = bim_.index().keys();
    bim_chrom_ = bim_["chrom"].as<std::string>();
    bim_pos_ = bim_["pos"].as<std::int32_t>();
    bim_a1_ = bim_["A1"].as<std::string>();
    bim_a2_ = bim_["A2"].as<std::string>();
}

SnpEffectsWriter::~SnpEffectsWriter() = default;

auto SnpEffectsWriter::write() -> void
{
    additive_ = result_->genetic(GeneticMode::A);
    dominant_ = result_->genetic(GeneticMode::D);
    if (additive_ == nullptr)
    {
        return;
    }

    write_header();

    for (Index i = 0; i < additive_->coeffs.size(); ++i)
    {
        row_buf_.clear();
        write_snp_row(i);
        writer_->write(row_buf_);
    }
}

auto SnpEffectsWriter::write_header() -> void
{
    Index n_add_components = 0;
    Index n_dom_components = 0;
    if (additive_ != nullptr && additive_->group)
    {
        n_add_components = assignment(*additive_->group).comp_probs.cols();
    }
    if (dominant_ != nullptr && dominant_->group)
    {
        n_dom_components = assignment(*dominant_->group).comp_probs.cols();
    }

    std::string h
        = "Index\tID\tChrom\tPosition\tA1\tA2\tA1Freq"
          "\tAdd\tAddSE\tAddPVE";

    if (n_add_components > 2)
    {
        for (Index comp = 0; comp < n_add_components; ++comp)
        {
            h += fmt::format("\tpi_{}", comp);
        }
    }
    h += "\tPIP";

    if (dominant_ != nullptr)
    {
        h += "\tDom\tDomSE\tDomPVE";
        if (n_dom_components > 2)
        {
            for (Index comp = 0; comp < n_dom_components; ++comp)
            {
                h += fmt::format("\tdpi_{}", comp);
            }
        }
        h += "\tdPIP";
    }

    writer_->write(h);
}

auto SnpEffectsWriter::write_snp_row(Index snp_index) -> void
{
    row_buf_ += fmt::format("{}\t", snp_index + 1);

    write_snp_basic_info(snp_index);
    write_effects(additive_, snp_index);
    write_component_probs(additive_, snp_index);
    write_pip(additive_, snp_index);
    write_effects(dominant_, snp_index);
    write_component_probs(dominant_, snp_index);
    write_pip(dominant_, snp_index);
}

auto SnpEffectsWriter::write_snp_basic_info(Index snp_index) -> void
{
    auto row = static_cast<std::size_t>(snp_index);
    if (snp_index < static_cast<Index>(bim_.rows()))
    {
        row_buf_ += fmt::format(
            "{}\t{}\t{}\t{}\t{}",
            bim_keys_[row],
            bim_chrom_[row],
            bim_pos_[row],
            bim_a1_[row],
            bim_a2_[row]);

        if (result_->allele_freq().size() > snp_index)
        {
            row_buf_
                += fmt::format("\t{:.6f}", result_->allele_freq()(snp_index));
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
    const GeneticSummary* effect,
    Index snp_index) -> void
{
    if (effect == nullptr || snp_index >= effect->coeffs.size())
    {
        return;
    }

    row_buf_ += fmt::format(
        "\t{:.6f}\t{:.6f}",
        effect->coeffs.mean(snp_index),
        effect->coeffs.stddev(snp_index));

    if (effect->pve.size() > snp_index)
    {
        row_buf_ += fmt::format("\t{:.6e}", effect->pve.mean(snp_index));
    }
    else
    {
        row_buf_ += "\t0.0";
    }
}

auto SnpEffectsWriter::write_component_probs(
    const GeneticSummary* effect,
    Index snp_index) -> void
{
    if (effect == nullptr || !effect->group)
    {
        return;
    }
    const auto& base = assignment(*effect->group);
    if (base.comp_probs.cols() <= 2 || snp_index >= base.comp_probs.rows())
    {
        return;
    }

    for (Index comp = 0; comp < base.comp_probs.cols(); ++comp)
    {
        row_buf_ += fmt::format("\t{:.6f}", base.comp_probs(snp_index, comp));
    }
}

auto SnpEffectsWriter::write_pip(const GeneticSummary* effect, Index snp_index)
    -> void
{
    if (effect == nullptr)
    {
        return;
    }

    if (effect->group && assignment(*effect->group).pip.size() > snp_index)
    {
        row_buf_ += fmt::format(
            "\t{:.6f}", assignment(*effect->group).pip(snp_index));
    }
    else
    {
        row_buf_ += "\t1.0";
    }
}

}  // namespace gelex
