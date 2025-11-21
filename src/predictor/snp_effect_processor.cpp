#include "snp_effect_processor.h"

#include <expected>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

#include "Eigen/Core"
#include "data/loader.h"
#include "data/parser.h"

namespace gelex
{

auto SnpEffectProcessor::create(const std::filesystem::path& snp_eff_path)
    -> std::expected<SnpEffectProcessor, Error>
{
    auto snp_infos = parse_snp_eff_file(snp_eff_path.string());
    if (!snp_infos)
    {
        return std::unexpected(snp_infos.error());
    }

    Eigen::VectorXd add(snp_infos->size());
    Eigen::VectorXd dom(snp_infos->size());
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(snp_infos->size());
         ++i)
    {
        add(i) = (*snp_infos)[i].add_effect;
        dom(i) = (*snp_infos)[i].dom_effect;
    }
    bool has_dominant_effect = dom.array().isFinite().any();
    SnpEffectProcessor processor;
    processor.snp_effects_ = std::move(*snp_infos);
    processor.add_ = std::move(add);
    processor.dom_ = std::move(dom);
    processor.has_dominant_effects_ = has_dominant_effect;
    return processor;
}

std::expected<std::vector<SnpEffect>, Error>
SnpEffectProcessor::parse_snp_eff_file(const std::string& snp_eff_file)
{
    auto file = detail::open_file<std::ifstream>(snp_eff_file, std::ios::in);
    if (!file)
    {
        return std::unexpected(file.error());
    }

    std::vector<SnpEffect> snp_infos;
    std::string line;
    if (!std::getline(*file, line))
    {
        return std::unexpected(
            Error{ErrorCode::InvalidData, "Empty .snp.eff file"});
    }

    auto header = detail::parse_string(line);
    const auto indices = assign_column_indices(header);

    // Parse data rows
    while (std::getline(*file, line))
    {
        auto row = detail::parse_string(line);
        auto snp_info = create_snp_info(row, indices);
        if (snp_info)
        {
            snp_infos.push_back(*snp_info);
        }
        else
        {
            return std::unexpected(snp_info.error());
        }
    }
    return snp_infos;
}

std::vector<double> SnpEffectProcessor::calculate_total_genetic_value(
    const std::vector<std::vector<int>>& genotypes,
    const std::vector<SnpEffect>& snp_infos)
{
    if (genotypes.empty() || snp_infos.empty())
    {
        return {};
    }

    // Validate input dimensions
    size_t n_individuals = genotypes[0].size();
    size_t n_snps = genotypes.size();

    if (n_snps != snp_infos.size())
    {
        // Mismatch between genotype data and SNP effects
        return {};
    }

    // Initialize result vector
    std::vector<double> total_values(n_individuals, 0.0);

    // Calculate total genetic value for each individual
    for (size_t snp_idx = 0; snp_idx < n_snps; ++snp_idx)
    {
        const auto& snp_info = snp_infos[snp_idx];
        const auto& snp_genotypes = genotypes[snp_idx];

        for (size_t ind_idx = 0; ind_idx < n_individuals; ++ind_idx)
        {
            double gevi = calculate_gevi(snp_genotypes[ind_idx], snp_info);
            total_values[ind_idx] += gevi;
        }
    }

    return total_values;
}

ColumnIndices SnpEffectProcessor::assign_column_indices(
    std::span<const std::string_view> header_columns)
{
    ColumnIndices indices;

    for (int i = 0; i < static_cast<int>(header_columns.size()); ++i)
    {
        const auto& column = header_columns[i];
        if (column == "ID")
        {
            indices.id = i;
        }
        else if (column == "A1")
        {
            indices.a1 = i;
        }
        else if (column == "A2")
        {
            indices.a2 = i;
        }
        else if (column == "A1Frq")
        {
            indices.a1frq = i;
        }
        else if (column == "Add")
        {
            indices.add = i;
        }
        else if (column == "Dom")
        {
            indices.dom = i;
        }
    }

    return indices;
}

auto SnpEffectProcessor::create_snp_info(
    std::span<const std::string_view> columns,
    const ColumnIndices& indices) -> std::expected<SnpEffect, Error>
{
    SnpEffect info;
    info.meta.id = columns[indices.id];
    info.meta.a1 = columns[indices.a1].empty() ? ' ' : columns[indices.a1][0];
    info.meta.a2 = columns[indices.a2].empty() ? ' ' : columns[indices.a2][0];
    auto frq = detail::try_parse_double(columns[indices.a1frq]);
    if (!frq)
    {
        return std::unexpected(frq.error());
    }
    info.p_freq = *frq;

    auto add_effect = detail::try_parse_double(columns[indices.add]);
    if (!add_effect)
    {
        return std::unexpected(add_effect.error());
    }
    info.add_effect = *add_effect;

    if (indices.dom != -1)
    {
        auto dom_effect = detail::try_parse_double(columns[indices.dom]);
        if (!dom_effect)
        {
            return std::unexpected(dom_effect.error());
        }
        info.dom_effect = *dom_effect;
    }
    return info;
}

}  // namespace gelex
