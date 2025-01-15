#include "chenx/data/bed_reader.h"
#include <omp.h>
#include <cstdint>
#include <iostream>
#include "armadillo"

/**
 * @brief Construct a new BedReader object
 */
BedReader::BedReader(
    const std::string& bed_file,
    size_t chunk_size,
    bool dom,
    int threads)
    : bed_file_{bed_file}, chunk_size_{chunk_size}, threads_{threads}
{
    std::string base_path = bed_file_.substr(0, bed_file_.size() - 4);
    std::string fam_file = base_path + ".fam";
    std::string bim_file = base_path + ".bim";

    n_individuals_ = parseFam(fam_file);
    snps_ = parseBim(bim_file);
    OpenBed();

    bytes_per_snp_ = (n_individuals_ + 3) / 4;  // add 3 for correct rounding
    if (threads_ > 0)
    {
        omp_set_num_threads(threads_);
    }
    else
    {
        omp_set_num_threads(omp_get_max_threads());
    }
    m_geno_map = dom ? genotypeMap_dom : genotypeMap_add;
}

/**
 * @brief Destroy the BedReader object
 */
BedReader::~BedReader()
{
    if (fin_.is_open())
    {
        fin_.close();
    }
}

//------------------------------------------------------------------------------
// Private utilities
//------------------------------------------------------------------------------
uint64_t BedReader::parseFam(const std::string& fam_file)
{
    std::ifstream fin(fam_file);
    if (!fin.is_open())
    {
        throw std::runtime_error(
            "Error: Cannot open .fam file [" + fam_file + "].");
    }

    uint64_t count = 0;
    std::string line;
    while (std::getline(fin, line))
    {
        if (!line.empty())
        {
            count++;
        }
    }
    fin.close();
    return count;
}

std::vector<SNP> BedReader::parseBim(const std::string& bim_file)
{
    std::ifstream fin(bim_file);
    if (!fin.is_open())
    {
        throw std::runtime_error(
            "Error: Cannot open .bim file [" + bim_file + "].");
    }

    std::vector<SNP> snps;
    std::string line;
    while (std::getline(fin, line))
    {
        if (line.empty())
        {
            continue;
        }
        std::istringstream iss(line);
        SNP snp;
        std::string dummy;
        iss >> snp.chromosome >> dummy >> dummy >> snp.position >> snp.allele1
            >> snp.allele2;
        snps.emplace_back(snp);
    }
    fin.close();
    return snps;
}

void BedReader::OpenBed()
{
    // Ensure file extension is .bed

    fin_.open(bed_file_, std::ios::binary);
    if (!fin_.is_open())
    {
        throw std::runtime_error(
            "Error: Cannot open BED file [" + bed_file_ + "].");
    }

    char buffer_header[3];
    fin_.read(buffer_header, 3);

    // Check for correct BED file magic numbers
    if (fin_.gcount() != 3)
    {
        throw std::runtime_error("Error: BED file header is incomplete.");
    }
    if (buffer_header[0] != 0x6C || buffer_header[1] != 0x1B
        || (buffer_header[2] != 0x01))
    {
        throw std::runtime_error(
            "Error: Invalid BED file format. Magic numbers do not match. "
            "Make "
            "sure using Plink1.9 to do the conversion.");
    }
}

//------------------------------------------------------------------------------
// Public: Iterable interface
//------------------------------------------------------------------------------
bool BedReader::HasNext() const
{
    // We can read until we've processed _snps.size() SNPs
    // Each chunk processes _chunk_size SNPs, so hasNext() is true
    // if the next chunk start index hasn't exceeded the total SNP count.
    return (current_chunk_index_ < snps_.size());
}

arma::dmat BedReader::GetNextChunk()
{
    if (!HasNext())
    {
        // Return empty if no more data
        std::cerr
            << "Warning: No more data to read from BED file. Reload if needed."
            << "\n";
        return {};
    }

    // Determine how many SNPs to read in this chunk
    size_t total_snps = snps_.size();
    size_t snps_remaining = total_snps - current_chunk_index_;
    size_t current_chunk_size = std::min(chunk_size_, snps_remaining);

    // Allocate buffer for raw bytes
    const uint64_t chunk_bytes = current_chunk_size * bytes_per_snp_;
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

    arma::dmat genotype_matrix(
        n_individuals_, current_chunk_size, arma::fill::zeros);

#pragma omp parallel for schedule(dynamic)
    for (uint64_t snp_idx = 0; snp_idx < current_chunk_size; ++snp_idx)
    {
        const size_t offset = snp_idx * bytes_per_snp_;
        for (uint64_t byte_idx = 0; byte_idx < bytes_per_snp_; ++byte_idx)
        {
            auto byte_val
                = static_cast<unsigned char>(buffer[offset + byte_idx]);
            for (int bit = 0; bit < 4; ++bit)
            {
                uint64_t ind = byte_idx * 4 + bit;
                if (ind >= n_individuals_)
                {
                    break;
                }
                const int genotype_code
                    = (byte_val >> (2 * bit)) & 0x03;  // encode byte to int
                genotype_matrix.at(ind, snp_idx) = m_geno_map[genotype_code];
            }
        }
    }
    current_chunk_index_ += current_chunk_size;
    return genotype_matrix;
}
