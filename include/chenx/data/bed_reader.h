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
    // Disable copy constructors/assignment operators
    BedReader(const BedReader&) = delete;
    BedReader(BedReader&&) noexcept = default;
    BedReader& operator=(const BedReader&) = delete;
    BedReader& operator=(BedReader&&) noexcept = default;

    explicit BedReader(std::string_view, size_t chunk_size = 10000);

    ~BedReader();

    bool HasNext() const;
    arma::dmat ReadChunk();
    // num snps
    uint64_t num_snps() const { return snps_.size(); }
    // snp id
    const std::vector<std::string>& snps() const noexcept { return snps_; }
    // num individuals
    uint64_t num_individuals() const noexcept { return individuals_.size(); }
    // individuals id
    const std::vector<std::string>& individuals() const noexcept
    {
        return individuals_;
    }

    uint64_t current_chunk_index() const noexcept
    {
        return current_chunk_index_;
    }

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

    static constexpr std::array<double, 4> genotype_map = {2.0, 1.0, 1.0, 0.0};
    // std::string_view is not accepted by std::ifstream
    static std::vector<std::string> parseFam(const std::string& fam_file);
    static std::vector<std::string> parseBim(const std::string& bim_file);

    arma::dmat Decode(const std::vector<char>& buffer, uint64_t chunk_size);
    void OpenBed();
};
}  // namespace chenx
