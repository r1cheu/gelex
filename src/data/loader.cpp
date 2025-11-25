#include "loader.h"

#include <algorithm>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <string_view>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>

#include <fmt/ranges.h>
#include "../src/data/parser.h"
#include "gelex/exception.h"
#include "gelex/logger.h"

namespace gelex::detail
{

// =============================================================================
// Helper Function: Safe Open
// =============================================================================
static std::ifstream open_stream(const std::filesystem::path& path)
{
    return detail::open_file<std::ifstream>(path, std::ios_base::in);
}

// =============================================================================
// PhenotypeLoader
// =============================================================================

PhenotypeLoader::PhenotypeLoader(
    const std::filesystem::path& path,
    int pheno_column,
    bool iid_only)
{
    auto file = open_stream(path);

    std::string line;
    if (!std::getline(file, line))
    {
        throw InvalidFileException("Empty phenotype file");
    }

    auto header = parse_header(line);
    if (pheno_column < 2 || static_cast<size_t>(pheno_column) >= header.size())
    {
        throw InvalidRangeException(
            std::format("Phenotype column {} out of range", pheno_column));
    }
    name_ = std::string(header[pheno_column]);

    int n_line = 0;
    data_.reserve(1024);

    while (std::getline(file, line))
    {
        n_line++;
        if (line.empty())
        {
            continue;
        }
        try
        {
            data_.emplace(
                parse_id(line, iid_only), parse_nth_double(line, pheno_column));
        }
        catch (const std::exception& e)
        {
            throw InvalidDataException(
                std::format("{} (line {})", e.what(), n_line + 1));
        }
    }

    gelex::logging::get()->info(
        "Loaded {} samples with phenotype '{}'.", data_.size(), name_);
}

Eigen::VectorXd PhenotypeLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
{
    Eigen::VectorXd result(id_map.size());
    result.setConstant(std::numeric_limits<double>::quiet_NaN());

    for (const auto& [id, value] : data_)
    {
        if (auto it = id_map.find(id); it != id_map.end())
        {
            result(it->second) = value;
        }
    }
    return result;
}

// =============================================================================
// QcovarLoader
// =============================================================================

QcovarLoader::QcovarLoader(const std::filesystem::path& path, bool iid_only)
{
    auto file = open_stream(path);

    std::string line;
    if (!std::getline(file, line))
    {
        throw InvalidFileException("Empty qcovar file");
    }

    auto header_view = parse_header(line);
    if (header_view.size() < 3)
    {
        throw InvalidRangeException("Qcovar file must have > 2 columns");
    }

    names_.reserve(header_view.size() - 2);
    for (size_t i = 2; i < header_view.size(); ++i)
    {
        names_.emplace_back(header_view[i]);
    }

    int n_line = 0;
    data_.reserve(1024);

    while (std::getline(file, line))
    {
        n_line++;
        if (line.empty())
        {
            continue;
        }

        try
        {
            data_.emplace(parse_id(line, iid_only), parse_all_doubles(line, 2));
        }
        catch (const std::exception& e)
        {
            throw InvalidDataException(
                std::format("{} (line {})", e.what(), n_line + 1));
        }
    }

    gelex::logging::get()->info(
        "Loaded {} samples with {} qcovars.", data_.size(), names_.size());
}

Eigen::MatrixXd QcovarLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
{
    const auto n_samples = static_cast<Eigen::Index>(id_map.size());
    const auto n_covars = static_cast<Eigen::Index>(names_.size());

    Eigen::MatrixXd result(n_samples, n_covars);
    result.setZero();

    for (const auto& [id, values] : data_)
    {
        if (auto it = id_map.find(id); it != id_map.end())
        {
            const Eigen::Index row_idx = it->second;
            const size_t copy_len
                = std::min(static_cast<size_t>(n_covars), values.size());
            if (copy_len > 0)
            {
                std::copy_n(
                    values.begin(), copy_len, result.row(row_idx).data());
            }
        }
    }
    return result;
}

// =============================================================================
// CovarLoader (Categorical)
// =============================================================================

CovarLoader::CovarLoader(const std::filesystem::path& path, bool iid_only)
{
    auto file = open_stream(path);

    std::string line;
    if (!std::getline(file, line))
    {
        throw InvalidFileException("Empty covar file");
    }

    auto header_view = parse_header(line);
    if (header_view.size() < 3)
    {
        throw InvalidRangeException("Covar file too few columns");
    }

    names_.reserve(header_view.size() - 2);
    for (size_t i = 2; i < header_view.size(); ++i)
    {
        names_.emplace_back(header_view[i]);
    }

    int n_line = 0;
    while (std::getline(file, line))
    {
        n_line++;
        if (line.empty())
        {
            continue;
        }
        try
        {
            raw_data_.emplace(
                parse_id(line, iid_only),
                std::ranges::to<std::vector<std::string>>(
                    parse_string(line, 2)));
        }
        catch (const std::exception& e)
        {
            throw InvalidDataException(
                std::format("{} (line {})", e.what(), n_line + 1));
        }
    }

    build_optimized_encodings();

    gelex::logging::get()->info(
        "Loaded {} samples with {} categorical covars.",
        raw_data_.size(),
        names_.size());
}

void CovarLoader::build_optimized_encodings()
{
    // 1. 收集每个协变量的所有唯一值 (Levels)
    std::vector<std::unordered_set<std::string_view>> levels_per_col(
        names_.size());

    for (const auto& [id, row_values] : raw_data_)
    {
        for (size_t i = 0; i < row_values.size() && i < levels_per_col.size();
             ++i)
        {
            if (!row_values[i].empty())
            {
                levels_per_col[i].insert(row_values[i]);
            }
        }
    }

    // 2. 构建编码映射
    optimized_encodings_.resize(names_.size());
    total_dummy_vars_ = 0;

    for (size_t i = 0; i < names_.size(); ++i)
    {
        const auto& levels = levels_per_col[i];
        if (levels.size() < 2)
        {
            continue;  // 只有一个 level 或空，不产生 dummy 变量
        }

        // 排序以保证确定性
        std::vector<std::string> sorted_levels(levels.begin(), levels.end());
        std::ranges::sort(sorted_levels);

        // Reference level method: drop first level
        // Level[0] -> [0, 0, ...], Level[1] -> [1, 0, ...], Level[2] -> [0, 1,
        // ...]
        size_t n_dummies = sorted_levels.size() - 1;

        // 第一个 level 全 0
        EncodedCovariate ref_enc;
        ref_enc.encoding.assign(n_dummies, 0);
        ref_enc.start_col_offset = static_cast<int>(total_dummy_vars_);
        optimized_encodings_[i].emplace(sorted_levels[0], ref_enc);

        for (size_t j = 1; j < sorted_levels.size(); ++j)
        {
            EncodedCovariate enc;
            enc.encoding.assign(n_dummies, 0);
            enc.encoding[j - 1] = 1;
            enc.start_col_offset = static_cast<int>(total_dummy_vars_);

            optimized_encodings_[i].emplace(sorted_levels[j], std::move(enc));
        }
        total_dummy_vars_ += n_dummies;
    }
}

Eigen::MatrixXd CovarLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
{
    Eigen::MatrixXd result(id_map.size(), total_dummy_vars_);
    result.setZero();

    for (const auto& [id, raw_values] : raw_data_)
    {
        auto id_it = id_map.find(id);
        if (id_it == id_map.end())
        {
            continue;
        }

        const Eigen::Index row_idx = id_it->second;

        const size_t n_cols_process
            = std::min(raw_values.size(), optimized_encodings_.size());

        for (size_t i = 0; i < n_cols_process; ++i)
        {
            // 查找当前值对应的编码信息
            const auto& encoding_map = optimized_encodings_[i];
            if (encoding_map.empty())
            {
                continue;
            }

            if (auto enc_it = encoding_map.find(raw_values[i]);
                enc_it != encoding_map.end())
            {
                const auto& info = enc_it->second;
                // 只有非零向量才需要写入 (优化稀疏写入)
                // 实际上 info.encoding 对于 ref level 全是 0，对于其他只有一个
                // 1 我们可以进一步优化：只存哪个位置是 1，不需要拷贝整个 vector

                // 通用写法：
                const Eigen::Index offset = info.start_col_offset;
                const auto len
                    = static_cast<Eigen::Index>(info.encoding.size());

                // 利用 Eigen 的 map 直接赋值，避免循环
                if (len > 0)
                {
                    // cast<double> 是必须的，因为 encoding 是 int
                    result.row(row_idx).segment(offset, len)
                        = Eigen::Map<const Eigen::VectorXi>(
                              info.encoding.data(), len)
                              .cast<double>();
                }
            }
        }
    }
    return result;
}

// =============================================================================
// BimLoader
// =============================================================================

BimLoader::BimLoader(const std::filesystem::path& path)
{
    auto file = open_stream(path);

    // estimate reserve size
    auto fsize = std::filesystem::file_size(path);
    snp_meta_.reserve(fsize / 30);
    char delimiter = detect_delimiter(file);

    std::string line;
    int n_line = 0;
    std::array<std::string_view, 6> cols;

    while (std::getline(file, line))
    {
        n_line++;
        if (line.empty())
        {
            continue;
        }

        auto tokens = line | std::views::split(delimiter);

        size_t col_idx = 0;
        for (auto&& sub_range : tokens)
        {
            if (sub_range.empty())
            {
                continue;
            }
            if (col_idx < 6)
            {
                cols[col_idx++] = std::string_view(sub_range);
            }
            else
            {
                break;
            }
        }

        if (col_idx < 6)
        {
            throw InconsistentColumnCountException(
                std::format("Bim line {} has fewer than 6 columns", n_line));
        }

        SnpMeta meta;
        meta.chrom = cols[0];
        meta.id = cols[1];

        auto res = std::from_chars(
            cols[3].data(), cols[3].data() + cols[3].size(), meta.position);

        if (res.ec != std::errc())
        {
            throw InvalidDataException(
                std::format("Invalid pos at line {}", n_line));
        }

        meta.A1 = cols[4][0];
        meta.A2 = cols[5][0];

        snp_meta_.push_back(std::move(meta));
    }
}

char BimLoader::detect_delimiter(std::ifstream& file)
{
    char delimiter = '\t';
    {
        std::string probe_line;
        while (std::getline(file, probe_line) && probe_line.empty())
        {
        }

        if (!probe_line.empty())
        {
            if (probe_line.find('\t') == std::string::npos)
            {
                delimiter = ' ';
            }
        }

        file.clear();
        file.seekg(0, std::ios::beg);
    }
    return delimiter;
}

std::vector<std::string> BimLoader::get_ids() const
{
    std::vector<std::string> ids;
    ids.reserve(snp_meta_.size());
    for (const auto& s : snp_meta_)
    {
        ids.push_back(s.id);
    }
    return ids;
}

// =============================================================================
// FamLoader
// =============================================================================

FamLoader::FamLoader(const std::filesystem::path& path, bool iid_only)
{
    auto file = open_stream(path);

    std::string line;
    int n_line = 0;
    ids_.reserve(1024);  // Start small

    while (std::getline(file, line))
    {
        n_line++;
        if (line.empty())
        {
            continue;
        }
        try
        {
            ids_.emplace_back(parse_id(line, iid_only, ' '));
        }
        catch (const std::exception& e)
        {
            throw InvalidDataException(
                std::format(
                    "{} (line {}) at {}", e.what(), n_line, path.string()));
        }
    }

    data_.reserve(ids_.size());
    for (size_t i = 0; i < ids_.size(); ++i)
    {
        data_[ids_[i]] = static_cast<Eigen::Index>(i);
    }

    gelex::logging::get()->info(
        "Loaded {} samples from fam file.", ids_.size());
}

}  // namespace gelex::detail
