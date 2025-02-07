#pragma once

#include <omp.h>
#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include <armadillo>

namespace chenx
{

// Read Bed file
// Example:
// BedReader reader("test.bed", 1000)
// while (reader.HasNext()) {
//     auto genotype_mat {reader.ReadChunk()}
// }
class BedReader
{
   public:
    // Disable copy and move constructors/assignment operators
    BedReader(const BedReader&) = delete;
    BedReader(BedReader&&) = delete;
    BedReader& operator=(const BedReader&) = delete;
    BedReader& operator=(BedReader&&) = delete;

    explicit BedReader(const std::string& bed_file, size_t chunk_size = 10000);

    ~BedReader();

    bool HasNext() const;
    arma::dmat ReadChunk();
    // num snps
    uint64_t n_snps() const { return snps_.size(); }
    // snp id
    const std::vector<std::string>& snps() const { return snps_; }
    // num individuals
    uint64_t n_individuals() const { return individuals_.size(); }
    // individuals id
    const std::vector<std::string>& individuals() const { return individuals_; }

    uint64_t current_chunk_index() const { return current_chunk_index_; }

   private:
    std::ifstream fin_;
    std::string bed_file_;
    std::string bim_file_;
    std::string fam_file_;

    std::vector<std::string> snps_;
    std::vector<std::string> individuals_;

    uint64_t chunk_size_;
    uint64_t current_chunk_index_{};
    uint64_t current_chunk_size_{};
    uint64_t bytes_per_snp_{};

    static constexpr std::array<double, 4> genotype_map = {0.0, 1.0, 1.0, 2.0};
    static std::vector<std::string> parseFam(const std::string& fam_file);
    static std::vector<std::string> parseBim(const std::string& bim_file);

    arma::dmat Decode(const std::vector<char>& buffer, uint64_t chunk_size);
    void OpenBed();
};
}  // namespace chenx
