#pragma once

#include <array>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

static constexpr Eigen::Index DEFAULT_CHUNK_SIZE = 10000;
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
        Eigen::Index chunk_size = DEFAULT_CHUNK_SIZE,
        const std::vector<std::string>& target_order = {});

    BedReader(const BedReader&) = delete;
    BedReader(BedReader&&) noexcept = default;
    BedReader& operator=(const BedReader&) = delete;
    BedReader& operator=(BedReader&&) noexcept = default;

    ~BedReader();

    Eigen::MatrixXd read_chunk();
    void reset();
    bool has_next() const;

    Eigen::Index num_individuals() const { return individuals_.size(); }
    Eigen::Index num_snps() const { return snps_.size(); }
    Eigen::Index chunk_size() const { return chunk_size_; }
    const std::vector<std::string>& snps() const { return snps_; }
    const std::vector<std::string>& individuals() const { return individuals_; }
    Eigen::Index current_chunk_index() const { return current_chunk_index_; }
    Eigen::Index current_chunk_size() const { return current_chunk_size_; }

   private:
    void open_bed();
    void seek_to_bed_start();

    void parse_fam(
        const std::string& fam_file,
        const std::vector<std::string>& target_order);

    std::vector<std::string> parse_bim(const std::string& bim_file);
    Eigen::MatrixXd decode(
        const std::vector<char>& buffer,
        Eigen::Index chunk_size);

    std::ifstream fin_;
    std::string bed_file_;
    std::string bim_file_;
    std::string fam_file_;

    std::vector<std::string> snps_;
    std::vector<std::string> individuals_;

    std::vector<Eigen::Index> file_index_to_target_index_;
    std::vector<bool> file_index_is_kept_;
    Eigen::Index total_samples_in_file_{};

    Eigen::Index chunk_size_;
    Eigen::Index current_chunk_index_{};
    Eigen::Index current_chunk_size_{};
    Eigen::Index bytes_per_snp_{};

    static constexpr std::array<double, 4> add_map = {2.0, 1.0, 1.0, 0.0};
};
}  // namespace gelex
