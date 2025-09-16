#include "gelex/data/bed_reader.h"

#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>

namespace gelex
{
using Eigen::Index;

char detect_separator(const std::string& line)
{
    Index tab_count = 0;
    Index space_count = 0;

    for (char c : line)
    {
        if (c == '\t')
        {
            tab_count++;
        }
        else if (c == ' ')
        {
            space_count++;
        }
    }

    // If we have tabs, assume tab-separated
    return (tab_count > 0) ? '\t' : ' ';
}

std::string find_second(const std::string& line, char separator)
{
    auto first = line.find(separator) + 1;
    auto second = line.find(separator, first);
    return line.substr(first, second - first);
}

BedReader::BedReader(
    std::string_view bed_file,
    Index chunk_size,
    const std::vector<std::string>& target_order)
    : bed_file_{bed_file}, chunk_size_{chunk_size}
{
    std::string base_path = bed_file_.substr(0, bed_file_.size() - 4);
    std::string fam_file = base_path + ".fam";
    std::string bim_file = base_path + ".bim";

    parse_fam(fam_file, target_order);
    snps_ = parse_bim(bim_file);

    open_bed();
    bytes_per_snp_
        = (total_samples_in_file_ + 3) / 4;  // add 3 for correct rounding
}

BedReader::~BedReader()
{
    if (fin_.is_open())
    {
        fin_.close();
    }
}

void BedReader::parse_fam(
    const std::string& fam_file,
    const std::vector<std::string>& target_order)
{
    std::ifstream fin(fam_file);
    if (!fin.is_open())
    {
        throw std::runtime_error(
            "Error: Cannot open .fam file [" + fam_file + "].");
    }
    std::vector<std::string> ids_in_file_order;
    std::string line;
    while (std::getline(fin, line))
    {
        if (!line.empty())
        {
            char separator = detect_separator(line);
            ids_in_file_order.emplace_back(find_second(line, separator));
        }
    }
    fin.close();

    total_samples_in_file_ = ids_in_file_order.size();
    if (total_samples_in_file_ == 0)
    {
        return;
    }

    file_index_is_kept_.assign(total_samples_in_file_, false);
    file_index_to_target_index_.assign(total_samples_in_file_, -1);

    if (target_order.empty())
    {
        individuals_ = ids_in_file_order;
        for (Index i = 0; i < total_samples_in_file_; ++i)
        {
            file_index_is_kept_[i] = true;
            file_index_to_target_index_[i] = i;
        }
    }
    else
    {
        individuals_ = target_order;
        std::unordered_map<std::string, Index> target_id_to_idx;
        for (Index i = 0; i < target_order.size(); ++i)
        {
            target_id_to_idx[target_order[i]] = i;
        }

        std::vector<std::string> found_targets;
        found_targets.reserve(target_order.size());

        for (Index i = 0; i < total_samples_in_file_; ++i)
        {
            const auto& current_id = ids_in_file_order[i];
            auto it = target_id_to_idx.find(current_id);

            if (it != target_id_to_idx.end())
            {
                auto target_idx = it->second;
                file_index_is_kept_[i] = true;
                file_index_to_target_index_[i] = target_idx;
                found_targets.push_back(current_id);
            }
        }

        if (found_targets.size() != target_order.size())
        {
            std::unordered_set<std::string> found_set(
                found_targets.begin(), found_targets.end());
            std::string missing_ids;
            for (const auto& id : target_order)
            {
                if (!found_set.contains(id))
                {
                    missing_ids += " " + id;
                }
            }
            throw std::runtime_error(
                "Error: The following target individuals were not found in the "
                ".fam file:"
                + missing_ids);
        }
    }
}

std::vector<std::string> BedReader::parse_bim(const std::string& bim_file)
{
    std::ifstream fin(bim_file);
    if (!fin.is_open())
    {
        throw std::runtime_error(
            "Error: Cannot open .bim file [" + bim_file + "].");
    }

    std::vector<std::string> snps;
    std::string line;

    while (std::getline(fin, line))
    {
        if (line.empty())
        {
            continue;
        }
        char separator = detect_separator(line);
        snps.emplace_back(find_second(line, separator));
    }
    fin.close();
    return snps;
}

void BedReader::open_bed()
{
    fin_.open(bed_file_, std::ios::binary);
    if (!fin_.is_open())
    {
        throw std::runtime_error(
            "Error: Cannot open BED file [" + bed_file_ + "].");
    }

    std::array<char, 3> buffer_header{};
    fin_.read(buffer_header.data(), 3);
    if (buffer_header[0] != 0x6C || buffer_header[1] != 0x1B
        || (buffer_header[2] != 0x01))
    {
        throw std::runtime_error(
            "Error: Invalid BED file format. Magic numbers do not match. "
            "Make "
            "sure using Plink1.9 to do the conversion.");
    }
}

bool BedReader::has_next() const
{
    return current_chunk_index() < num_snps();
}

void BedReader::reset()
{
    if (!fin_.is_open())
    {
        open_bed();
    }
    else
    {
        // Clear any error flags
        fin_.clear();

        // Return to the start of data (after the header)
        seek_to_bed_start();

        // reset the tracking variables
        current_chunk_index_ = 0;
        current_chunk_size_ = 0;
    }
}

void BedReader::seek_to_bed_start()
{
    // Seek to the beginning of file + 3 bytes (skip the magic number)
    fin_.seekg(3, std::ios::beg);
}

Eigen::MatrixXd BedReader::read_chunk()
{
    if (!has_next())
    {
        std::cerr
            << "Warning: No more data to read from BED file. Reload if needed."
            << "\n";
        return Eigen::MatrixXd();
    }

    current_chunk_size_
        = std::min(chunk_size_, num_snps() - current_chunk_index());

    const Index chunk_bytes = current_chunk_size_ * bytes_per_snp_;
    std::vector<char> buffer(chunk_bytes);

    fin_.read(buffer.data(), static_cast<int64_t>(chunk_bytes));
    Index bytes_read = fin_.gcount();
    if (bytes_read != chunk_bytes)
    {
        throw std::runtime_error(
            "Error: Incomplete SNP data in BED file. "
            "Expected "
            + std::to_string(chunk_bytes) + " bytes, got "
            + std::to_string(bytes_read));
    }

    Eigen::MatrixXd genotype_matrix = decode(buffer, current_chunk_size_);
    current_chunk_index_ += current_chunk_size_;
    return genotype_matrix;
}
Eigen::MatrixXd BedReader::decode(
    const std::vector<char>& buffer,
    Index chunk_size)
{
    Eigen::MatrixXd genotype_matrix
        = Eigen::MatrixXd::Zero(num_individuals(), chunk_size);

    for (Index snp_idx = 0; snp_idx < chunk_size; ++snp_idx)
    {
        const Index offset = snp_idx * bytes_per_snp_;
        for (Index byte_idx = 0; byte_idx < bytes_per_snp_; ++byte_idx)
        {
            auto byte_val
                = static_cast<unsigned char>(buffer[offset + byte_idx]);
            for (unsigned int bit = 0; bit < 4; ++bit)
            {
                Index file_index = (byte_idx * 4) + bit;

                if (file_index >= total_samples_in_file_)
                {
                    break;
                }

                if (file_index_is_kept_[file_index])
                {
                    Index target_row = file_index_to_target_index_[file_index];

                    unsigned int genotype_code
                        = (byte_val >> (2U * bit)) & 0x03U;
                    genotype_matrix(target_row, snp_idx)
                        = add_map[genotype_code];
                }
            }
        }
    }
    return genotype_matrix;
}

}  // namespace gelex
