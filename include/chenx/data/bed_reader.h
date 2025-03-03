#pragma once

#include <omp.h>
#include <array>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

#include <armadillo>

namespace chenx
{

static constexpr uint64_t DEFAULT_CHUNK_SIZE = 10000;
// Read Bed file
// Example:
// BedReader reader("test.bed", 1000)
// while (reader.HasNext()) {
//     auto genotype_mat {reader.ReadChunk()}
// }
//

std::string find_second(std::string& snps_line);

using strings = std::vector<std::string>;
class BedReader
{
   public:
    explicit BedReader(
        std::string_view,
        size_t chunk_size = DEFAULT_CHUNK_SIZE,
        const strings& dropped_individuals = {});

    BedReader(const BedReader&) = delete;
    BedReader(BedReader&&) noexcept = default;
    BedReader& operator=(const BedReader&) = delete;
    BedReader& operator=(BedReader&&) noexcept = default;

    ~BedReader();

    void Reset();

    uint64_t chunk_size() const noexcept { return chunk_size_; }
    bool HasNext() const;
    arma::dmat ReadChunk();
    uint64_t num_snps() const { return snps_.size(); }
    const strings& snps() const noexcept { return snps_; }
    uint64_t num_individuals() const noexcept { return individuals_.size(); }
    const strings& individuals() const noexcept { return individuals_; }
    uint64_t current_chunk_index() const noexcept
    {
        return current_chunk_index_;
    }
    uint64_t current_chunk_size() const noexcept { return current_chunk_size_; }

   private:
    std::ifstream fin_;
    std::string bed_file_;
    std::string bim_file_;
    std::string fam_file_;

    strings snps_;
    strings individuals_;

    std::unordered_set<uint64_t> exclude_index_;

    uint64_t chunk_size_;
    uint64_t current_chunk_index_{};
    uint64_t current_chunk_size_{};
    uint64_t bytes_per_snp_{};

    static constexpr std::array<double, 4> genotype_map = {2.0, 1.0, 1.0, 0.0};
    // std::string_view is not accepted by std::ifstream
    std::vector<std::string> parseFam(
        const std::string& fam_file,
        const strings& dropped_individuals);
    static std::vector<std::string> parseBim(const std::string& bim_file);

    arma::dmat Decode(const std::vector<char>& buffer, uint64_t chunk_size);
    void OpenBed();
    void SeekToBedStart();
};
}  // namespace chenx
