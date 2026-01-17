#include "grm_loader.h"
#include <fmt/format.h>

#include <algorithm>
#include <format>
#include <fstream>
#include <string>
#include <system_error>

#include "gelex/exception.h"
#include "parser.h"
#include "types/freq_effect.h"

namespace
{
auto get_type(std::string_view grm_path_stem) -> gelex::freq::GrmType
{
    if (grm_path_stem.contains("add"))
    {
        return gelex::freq::GrmType::A;
    }
    if (grm_path_stem.contains("dom"))
    {
        return gelex::freq::GrmType::D;
    }
    return gelex::freq::GrmType::Unknown;
}
}  // namespace
namespace gelex::detail
{

GrmLoader::GrmLoader(const std::filesystem::path& prefix)
    : bin_path_(prefix.string() + ".grm.bin"),
      id_path_(prefix.string() + ".grm.id"),
      type_(get_type(prefix.string()))
{
    load_sample_ids();
    init_mmap();
}

auto GrmLoader::load_sample_ids() -> void
{
    auto file = open_file<std::ifstream>(id_path_, std::ios::in);

    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }

        // Format: FID\tIID -> convert to FID_IID
        auto tab_pos = line.find('\t');
        if (tab_pos == std::string::npos)
        {
            // No tab found, use the line as both FID and IID
            sample_ids_.push_back(fmt::format("{}_{}", line, line));
        }
        else
        {
            auto fid = line.substr(0, tab_pos);
            auto iid = line.substr(tab_pos + 1);
            sample_ids_.push_back(std::format("{}_{}", fid, iid));
        }
    }

    if (sample_ids_.empty())
    {
        throw FileFormatException(
            std::format("{}: no sample IDs found", id_path_.string()));
    }

    num_samples_ = static_cast<Eigen::Index>(sample_ids_.size());
}

auto GrmLoader::init_mmap() -> void
{
    std::error_code ec;
    mmap_.map(bin_path_.string(), ec);
    if (ec)
    {
        throw FileOpenException(
            std::format("{}: failed to mmap file", bin_path_.string()));
    }

    // GRM binary format: lower triangle stored as float32
    // n * (n + 1) / 2 elements
    size_t expected_elements = static_cast<size_t>(num_samples_)
                               * (static_cast<size_t>(num_samples_) + 1) / 2;
    size_t expected_size = expected_elements * sizeof(float);

    if (mmap_.size() != expected_size)
    {
        throw FileFormatException(
            std::format(
                "{}: file size mismatch. Expected {} bytes ({} samples), got "
                "{} bytes",
                bin_path_.string(),
                expected_size,
                num_samples_,
                mmap_.size()));
    }
}

auto GrmLoader::load() const -> Eigen::MatrixXd
{
    Eigen::MatrixXd grm(num_samples_, num_samples_);

    const auto* data = reinterpret_cast<const float*>(mmap_.data());

    // Fill lower triangle and mirror to upper triangle
    for (Eigen::Index i = 0; i < num_samples_; ++i)
    {
        for (Eigen::Index j = 0; j <= i; ++j)
        {
            size_t idx = lower_triangle_index(i, j);
            auto value = static_cast<double>(data[idx]);
            grm(i, j) = value;
            grm(j, i) = value;
        }
    }

    return grm;
}

auto GrmLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
    -> Eigen::MatrixXd
{
    if (id_map.empty())
    {
        return Eigen::MatrixXd();
    }

    // Build file_id -> file_index mapping
    std::unordered_map<std::string, Eigen::Index> file_id_to_idx;
    file_id_to_idx.reserve(sample_ids_.size());
    for (Eigen::Index i = 0; i < num_samples_; ++i)
    {
        file_id_to_idx[sample_ids_[i]] = i;
    }

    // Build source_idx -> target_idx mapping
    // and find max target index to determine output matrix size
    Eigen::Index max_target_idx = 0;
    std::vector<std::pair<Eigen::Index, Eigen::Index>> idx_mapping;
    idx_mapping.reserve(id_map.size());

    for (const auto& [id, target_idx] : id_map)
    {
        auto it = file_id_to_idx.find(id);
        if (it == file_id_to_idx.end())
        {
            throw InvalidInputException(
                std::format(
                    "{}: sample ID '{}' not found in GRM file",
                    bin_path_.string(),
                    id));
        }

        idx_mapping.emplace_back(it->second, target_idx);
        max_target_idx = std::max(max_target_idx, target_idx);
    }

    // Allocate output matrix
    Eigen::Index out_size = max_target_idx + 1;
    Eigen::MatrixXd grm = Eigen::MatrixXd::Zero(out_size, out_size);

    const auto* data = reinterpret_cast<const float*>(mmap_.data());

    // Fill the matrix using index mapping
    for (size_t ii = 0; ii < idx_mapping.size(); ++ii)
    {
        auto [src_i, tgt_i] = idx_mapping[ii];

        for (size_t jj = 0; jj < idx_mapping.size(); ++jj)
        {
            auto [src_j, tgt_j] = idx_mapping[jj];

            // Read from lower triangle (ensure src_i >= src_j)
            Eigen::Index file_i = src_i;
            Eigen::Index file_j = src_j;
            if (file_i < file_j)
            {
                std::swap(file_i, file_j);
            }

            size_t idx = lower_triangle_index(file_i, file_j);
            grm(tgt_i, tgt_j) = static_cast<double>(data[idx]);
        }
    }

    return grm;
}

}  // namespace gelex::detail
