#pragma once

#include <array>
#include <string>
#include <unordered_set>
#include <vector>

#include <armadillo>

namespace gelex
{

static constexpr size_t DEFAULT_CHUNK_SIZE = 10000;
// Read Bed file
// Example:
// BedReader reader("test.bed", 1000)
// while (reader.has_next()) {
//     auto genotype_mat {reader.read_chunk()}
// }
//

std::string find_second(std::string& snps_line);

class BedReader
{
   public:
    explicit BedReader(
        std::string_view,
        size_t chunk_size = DEFAULT_CHUNK_SIZE,
        const std::vector<std::string>& dropped_ids = {});

    BedReader(const BedReader&) = delete;
    BedReader(BedReader&&) noexcept = default;
    BedReader& operator=(const BedReader&) = delete;
    BedReader& operator=(BedReader&&) noexcept = default;

    ~BedReader();

    void reset();

    size_t chunk_size() const noexcept { return chunk_size_; }
    bool has_next() const;
    arma::dmat read_chunk(bool add = true);
    size_t num_snps() const { return snps_.size(); }
    const std::vector<std::string>& snps() const noexcept { return snps_; }
    size_t num_individuals() const noexcept { return individuals_.size(); }
    const std::vector<std::string>& individuals() const noexcept
    {
        return individuals_;
    }
    size_t current_chunk_index() const noexcept { return current_chunk_index_; }
    size_t current_chunk_size() const noexcept { return current_chunk_size_; }

   private:
    std::ifstream fin_;
    std::string bed_file_;
    std::string bim_file_;
    std::string fam_file_;

    std::vector<std::string> snps_;
    std::vector<std::string> individuals_;

    std::unordered_set<size_t> exclude_index_;

    size_t chunk_size_;
    size_t current_chunk_index_{};
    size_t current_chunk_size_{};
    size_t bytes_per_snp_{};

    static constexpr std::array<double, 4> add_map = {2.0, 1.0, 1.0, 0.0};
    static constexpr std::array<double, 4> dom_map = {0.0, 1.0, 1.0, 0.0};

    std::vector<std::string> parse_fam(
        const std::string& fam_file,
        const std::vector<std::string>& dropped_ids);
    static std::vector<std::string> parse_bim(const std::string& bim_file);

    arma::dmat
    decode(const std::vector<char>& buffer, size_t chunk_size, bool add);
    void open_bed();
    void seek_to_bed_start();
};
}  // namespace gelex
