#include "chenx/data/bed_reader.h"

#include <omp.h>
#include <array>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include <armadillo>

namespace chenx
{
using strings = std::vector<std::string>;
BedReader::BedReader(
    std::string_view bed_file,
    strings&& dropped_individuals,
    size_t chunk_size)
    : bed_file_{bed_file},
      dropped_individuals_{std::move(dropped_individuals)},
      chunk_size_{chunk_size}
{
    std::string base_path = bed_file_.substr(0, bed_file_.size() - 4);
    std::string fam_file = base_path + ".fam";
    std::string bim_file = base_path + ".bim";

    individuals_ = parseFam(fam_file, dropped_individuals_);
    snps_ = parseBim(bim_file);

    OpenBed();
    bytes_per_snp_ = (num_individuals() + exclude_index_.size() + 3)
                     / 4;  // add 3 for correct rounding
}

BedReader::BedReader(
    std::string_view bed_file,
    const strings& dropped_individuals,
    size_t chunk_size)
    : BedReader{bed_file, strings{dropped_individuals}, chunk_size} {};

BedReader::~BedReader()
{
    if (fin_.is_open())
    {
        fin_.close();
    }
}

std::vector<std::string> BedReader::parseFam(
    const std::string& fam_file,
    const strings& dropped_individuals)
{
    std::unordered_set<std::string> exclude_set(
        dropped_individuals.begin(), dropped_individuals.end());

    std::ifstream fin(fam_file);
    if (!fin.is_open())
    {
        throw std::runtime_error(
            "Error: Cannot open .fam file [" + fam_file + "].");
    }

    std::string line;
    std::vector<std::string> individuals;
    uint64_t index{};
    while (std::getline(fin, line))
    {
        if (!line.empty())
        {
            auto first = line.find('\t') + 1;
            auto second = line.find('\t', first);
            auto individual = line.substr(first, second - first);
            if (!exclude_set.empty()
                && exclude_set.find(individual) != exclude_set.end())
            {
                exclude_index_.insert(index);
                ++index;
                continue;
            }
            individuals.emplace_back(individual);
            ++index;
        }
    }
    fin.close();
    return individuals;
}

std::vector<std::string> BedReader::parseBim(const std::string& bim_file)
{
    std::ifstream fin(bim_file);
    if (!fin.is_open())
    {
        throw std::runtime_error(
            "Error: Cannot open .bim file [" + bim_file + "].");
    }

    std::vector<std::string> snps;
    std::string line;
    std::string temp;

    while (std::getline(fin, line))
    {
        if (line.empty())
        {
            continue;
        }
        std::istringstream stream{line};
        std::string snp;
        for (int i{}; std::getline(stream, temp, '\t'); ++i)
        {
            switch (i)
            {
                case 0:
                    snp += temp;
                    break;
                case 1:
                case 2:
                    break;
                case 3:
                case 4:
                case 5:
                    snp += ":" + temp;
                    break;
                default:
                    break;
            }
        }
        snps.push_back(snp);
    }
    fin.close();
    return snps;
}

void BedReader::OpenBed()
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

bool BedReader::HasNext() const
{
    return current_chunk_index() < num_snps();
}

arma::dmat BedReader::ReadChunk()
{
    if (!HasNext())
    {
        std::cerr
            << "Warning: No more data to read from BED file. Reload if needed."
            << "\n";
        return {};
    }

    current_chunk_size_ = std::min(
        chunk_size_,
        num_snps() - current_chunk_index());  // chunk_size or snp remaining.

    const uint64_t chunk_bytes = current_chunk_size_ * bytes_per_snp_;
    std::vector<char> buffer(chunk_bytes);

    // Read from the current file position
    fin_.read(buffer.data(), static_cast<int64_t>(chunk_bytes));
    size_t bytes_read = fin_.gcount();
    if (bytes_read != chunk_bytes)
    {
        throw std::runtime_error(
            "Error: Incomplete SNP data in BED file. "
            "Expected "
            + std::to_string(chunk_bytes) + " bytes, got "
            + std::to_string(bytes_read));
    }

    arma::dmat genotype_matrix{Decode(buffer, current_chunk_size_)};
    current_chunk_index_ += current_chunk_size_;
    return genotype_matrix;
}

arma::dmat BedReader::Decode(
    const std::vector<char>& buffer,
    uint64_t chunk_size)
{
    arma::dmat genotype_matrix(
        num_individuals(), chunk_size, arma::fill::zeros);
    uint64_t num_individuals
        = individuals_.size()
          + exclude_index_.size();  // add back excluded individuals, since we
                                    // are read .bed

#pragma omp parallel for schedule(dynamic)
    for (uint64_t snp_idx = 0; snp_idx < chunk_size; ++snp_idx)
    {
        const size_t offset = snp_idx * bytes_per_snp_;
        for (uint64_t byte_idx = 0; byte_idx < bytes_per_snp_; ++byte_idx)
        {
            uint64_t adjust_individul_idx{};  // adjust for excluded individuals
            auto byte_val
                = static_cast<unsigned char>(buffer[offset + byte_idx]);
            for (unsigned int bit = 0; bit < 4; ++bit)
            {
                uint64_t ind = (byte_idx * 4) + bit;
                if (exclude_index_.find(ind) != exclude_index_.end())
                {
                    adjust_individul_idx++;
                    continue;
                }
                if (ind >= num_individuals)
                {
                    break;
                }
                unsigned int genotype_code
                    = (byte_val >> (2U * bit)) & 0x03U;  // encode byte to int
                genotype_matrix.at(ind - adjust_individul_idx, snp_idx)
                    = genotype_map[genotype_code];
            }
        }
    }
    return genotype_matrix;
}
}  // namespace chenx
