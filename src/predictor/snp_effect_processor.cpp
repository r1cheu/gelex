#include "snp_effect_processor.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace gelex
{

std::expected<std::vector<SnpInfo>, Error> SnpEffectProcessor::create(
    const std::string& snp_eff_file)
{
    return parse_snp_eff_file(snp_eff_file);
}

std::expected<std::vector<SnpInfo>, Error>
SnpEffectProcessor::parse_snp_eff_file(const std::string& snp_eff_file)
{
    std::ifstream file(snp_eff_file);
    if (!file.is_open())
    {
        return std::unexpected(
            Error{
                ErrorCode::FileNotFound,
                "Failed to open .snp.eff file: " + snp_eff_file});
    }

    std::vector<SnpInfo> snp_infos;
    std::string line;

    // Read and parse header
    if (!std::getline(file, line))
    {
        return std::unexpected(
            Error{ErrorCode::InvalidData, "Empty .snp.eff file"});
    }

    auto indices_result = parse_header(line);
    if (!indices_result)
    {
        return std::unexpected(indices_result.error());
    }
    const auto& indices = *indices_result;

    // Parse data rows
    while (std::getline(file, line))
    {
        auto snp_info = parse_snp_row(line, indices);
        if (snp_info)
        {
            snp_infos.push_back(*snp_info);
        }
    }

    if (snp_infos.empty())
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidData,
                "No valid SNP effects found in .snp.eff file"});
    }

    return snp_infos;
}

std::vector<double> SnpEffectProcessor::calculate_total_genetic_value(
    const std::vector<std::vector<int>>& genotypes,
    const std::vector<SnpInfo>& snp_infos)
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
    const std::vector<std::string>& header_columns)
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

std::expected<ColumnIndices, Error> SnpEffectProcessor::parse_header(
    const std::string& header_line)
{
    if (header_line.empty())
    {
        return std::unexpected(
            Error{ErrorCode::InvalidData, "Empty header line"});
    }

    // Parse header to understand column structure
    std::istringstream header_stream(header_line);
    std::vector<std::string> header_columns;
    std::string column;

    while (std::getline(header_stream, column, '\t'))
    {
        header_columns.push_back(column);
    }

    // Assign column indices using helper function
    ColumnIndices indices = assign_column_indices(header_columns);

    // Validate required columns
    if (!indices.has_required_columns())
    {
        return std::unexpected(
            Error{
                ErrorCode::WrongHeader,
                "Missing required columns in .snp.eff file header"});
    }

    return indices;
}

std::optional<SnpInfo> SnpEffectProcessor::parse_snp_row(
    const std::string& line,
    const ColumnIndices& indices)
{
    if (line.empty())
    {
        return std::nullopt;
    }

    std::istringstream line_stream(line);
    std::vector<std::string> columns;
    std::string value;

    while (std::getline(line_stream, value, '\t'))
    {
        columns.push_back(value);
    }

    return create_snp_info(columns, indices);
}

std::optional<SnpInfo> SnpEffectProcessor::create_snp_info(
    const std::vector<std::string>& columns,
    const ColumnIndices& indices)
{
    // Skip if we don't have enough columns
    if (columns.size() <= static_cast<size_t>(std::max(
            {indices.id, indices.a1, indices.a2, indices.a1frq, indices.add})))
    {
        return std::nullopt;
    }

    try
    {
        SnpInfo info;
        info.id = columns[indices.id];
        info.a1 = columns[indices.a1].empty() ? ' ' : columns[indices.a1][0];
        info.a2 = columns[indices.a2].empty() ? ' ' : columns[indices.a2][0];

        // Parse frequency
        if (columns[indices.a1frq] != "NA" && !columns[indices.a1frq].empty())
        {
            info.p_freq = std::stod(columns[indices.a1frq]);
        }
        else
        {
            // Skip SNPs with missing frequency
            return std::nullopt;
        }

        // Parse additive effect
        if (columns[indices.add] != "NA" && !columns[indices.add].empty())
        {
            info.add_effect = std::stod(columns[indices.add]);
        }
        else
        {
            // Skip SNPs with missing additive effect
            return std::nullopt;
        }

        // Parse dominant effect (optional)
        if (indices.dom != -1
            && columns.size() > static_cast<size_t>(indices.dom)
            && columns[indices.dom] != "NA" && !columns[indices.dom].empty())
        {
            info.dom_effect = std::stod(columns[indices.dom]);
        }
        else
        {
            info.dom_effect
                = 0.0;  // Default to zero if dominant effect not available
        }

        // Validate frequency range
        if (info.p_freq < 0.0 || info.p_freq > 1.0)
        {
            return std::nullopt;  // Skip invalid frequencies
        }

        return info;
    }
    catch (const std::exception& e)
    {
        // Skip malformed rows but continue processing
        return std::nullopt;
    }
}

}  // namespace gelex
