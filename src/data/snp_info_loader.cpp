#include "snp_info_loader.h"

#include <charconv>
#include <expected>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "../src/data/loader.h"
#include "../src/data/parser.h"

namespace gelex
{

auto SnpInfoLoader::create(const std::filesystem::path& bim_file_path)
    -> std::expected<SnpInfoLoader, Error>
{
    auto snp_info = read_bim_file(bim_file_path);
    if (!snp_info)
    {
        return std::unexpected(snp_info.error());
    }

    return SnpInfoLoader(std::move(*snp_info));
}

void SnpInfoLoader::set_allele_frequencies(const Eigen::VectorXd& frequencies)
{
    if (frequencies.size() != snp_info_.size())
    {
        throw std::invalid_argument(
            "Number of frequencies does not match number of SNPs");
    }

    for (size_t i = 0; i < snp_info_.size(); ++i)
    {
        snp_info_[i].a1_frq = frequencies(i);
    }
}

auto SnpInfoLoader::read_bim_file(const std::filesystem::path& path)
    -> std::expected<std::vector<SnpInfo>, Error>
{
    auto file = detail::open_file<std::ifstream>(path, std::ios_base::in);
    if (!file)
    {
        return std::unexpected(file.error());
    }

    constexpr static size_t bim_n_cols = 6;
    std::vector<SnpInfo> snp_info;
    std::string line;
    int n_line{};

    while (std::getline(*file, line))
    {
        if (line.empty())
        {
            continue;
        }

        if (detail::count_num_columns(line) != bim_n_cols)
        {
            return std::unexpected(
                Error{
                    ErrorCode::InconsistColumnCount,
                    std::format(
                        "BIM file must have exactly 6 columns (line {})",
                        n_line + 1)});
        }

        auto tokens = detail::parse_string(line, 0, "\t");
        if (tokens.size() < bim_n_cols)
        {
            return std::unexpected(
                Error{
                    ErrorCode::InvalidFile,
                    std::format(
                        "Failed to parse BIM file (line {})", n_line + 1)});
        }

        // Parse position using std::from_chars
        const std::string_view pos_str = tokens[3];
        int position = 0;
        auto result = std::from_chars(
            pos_str.data(), pos_str.data() + pos_str.size(), position);

        if (result.ec != std::errc() || position < 0)
        {
            return std::unexpected(
                Error{
                    ErrorCode::InvalidData,
                    std::format(
                        "Invalid position '{}' in BIM file (line {})",
                        pos_str,
                        n_line + 1)});
        }

        SnpInfo info;
        info.chrom = tokens[0];
        info.id = tokens[1];
        info.position = position;
        info.a1 = tokens[4];
        info.a2 = tokens[5];

        snp_info.push_back(std::move(info));
        n_line++;
    }

    return snp_info;
}

}  // namespace gelex
