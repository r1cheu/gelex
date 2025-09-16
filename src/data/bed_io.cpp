#include "gelex/data/bed_io.h"

#include <cassert>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <span>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <gelex/barkeep.h>
#include "gelex/data/io.h"

namespace gelex
{

using Eigen::Index;

using Eigen::VectorXd;

using index_map = std::unordered_map<std::string, Index>;

BedIO::BedIO(const std::string& bfile, bool iid_only)
    : bed_file_(bfile + ".bed"),
      fam_file_(bfile + ".fam"),
      bim_file_(bfile + ".bim")
{
    fam_ids_ = read_fam(fam_file_, iid_only);
    snp_ids_ = read_bim(bim_file_);
    create_map();
};

std::vector<std::string> BedIO::read_fam(
    const std::string& fam_path,
    bool iid_only)
{
    auto file = detail::open_or_throw<std::ifstream>(fam_path);
    std::string line;
    std::vector<std::string> fam_ids;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string fid;
        std::string iid;
        iss >> fid >> iid;
        if (iid_only)
        {
            fam_ids.emplace_back(iid);
        }
        else
        {
            fam_ids.emplace_back(fid.append("_").append(iid));
        }
    }
    return fam_ids;
}

std::vector<std::string> BedIO::read_bim(const std::string& bim_path)

{
    auto file = detail::open_or_throw<std::ifstream>(bim_path);
    std::string line;

    std::vector<std::string> snp_ids;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string snp_id;
        iss >> snp_id;
        iss >> snp_id;
        snp_ids.emplace_back(snp_id);

        std::string skip;
        while (iss >> skip)
        {
        }
    }
    return snp_ids;
}

std::ifstream BedIO::create_bed()
{
    auto bed = detail::open_or_throw<std::ifstream>(
        bed_file_, detail::file_type::binary);
    constexpr std::array<std::byte, 3> expected_bytes
        = {std::byte{0x6C}, std::byte{0x1B}, std::byte{0x01}};
    std::array<std::byte, 3> buffer{};
    bed.read(reinterpret_cast<char*>(buffer.data()), buffer.size());
    if (!bed)
    {
        throw std::runtime_error("Failed to read .bed header");
    }
    if (buffer != expected_bytes)
    {
        throw std::runtime_error("Invalid .bed header");
    }
    return bed;
}

void BedIO::read_locus(
    std::span<const std::byte> buffer,
    std::span<double> result)
{
    size_t write_count = 0;
    const size_t max_count = result.size();

    for (std::size_t i = 0; i < buffer.size(); ++i)
    {
        const auto value = std::to_integer<unsigned char>(buffer[i]);
        for (int j = 0; j < 4; ++j)
        {
            if (write_count >= max_count)
            {
                break;
            }

            int genotype = (value >> (j * 2)) & 0b11;
            double gen_value = genotype_map[genotype];
            result[write_count] = gen_value;

            write_count++;
        }
        if (write_count >= max_count)
        {
            break;
        }
    }
}

Eigen::VectorXd BedIO::rearange_locus(
    std::span<const Index> id_indices,
    std::span<const double> genotype)
{
    Eigen::VectorXd ordered_genotype(id_indices.size());
    for (Index i = 0; i < static_cast<Index>(id_indices.size()); ++i)
    {
        ordered_genotype(i) = genotype[id_indices[i]];
    }
    return ordered_genotype;
}

void BedIO::create_map()
{
    auto length = static_cast<Index>(fam_ids_.size());
    fam_map_.clear();
    fam_map_.reserve(fam_ids_.size());
    for (Index j = 0; j < length; ++j)
    {
        fam_map_.emplace(fam_ids_[j], j);
    }
}

std::vector<Index> BedIO::create_index_vector(
    std::span<const std::string> id_list)
{
    std::vector<Index> indices;
    indices.reserve(id_list.size());
    for (const auto& id : id_list)
    {
        auto it = fam_map_.find(id);
        if (it == fam_map_.end())
        {
            throw std::runtime_error("individual ID not found: " + id);
        }
        indices.emplace_back(it->second);
    }
    return indices;
}

}  // namespace gelex
