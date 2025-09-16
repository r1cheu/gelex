#pragma once

#include <cstddef>
#include <fstream>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

class BedIO
{
   public:
    explicit BedIO(const std::string& bfile, bool iid_only = false);

    size_t n_snp() const { return snp_ids_.size(); };
    size_t n_individuals() const { return fam_ids_.size(); }

    const auto& snp_ids() const { return snp_ids_; }
    const auto& fam_ids() const { return fam_ids_; }

    static std::vector<std::string> read_fam(
        const std::string& fam_path,
        bool iid_only);
    static std::vector<std::string> read_bim(const std::string& bim_path);

    std::ifstream create_bed();

    static void read_locus(
        std::span<const std::byte> buffer,
        std::span<double> result);

    static Eigen::VectorXd rearange_locus(
        std::span<const Eigen::Index> id_indices,
        std::span<const double> genotype);

    std::vector<Eigen::Index> create_index_vector(
        std::span<const std::string> id_list);

   private:
    static constexpr std::array<double, 4> genotype_map = {2, 1, 1, 0};
    void create_map();

    std::string bed_file_;
    std::string fam_file_;
    std::string bim_file_;

    std::unordered_map<std::string, Eigen::Index> fam_map_;

    std::vector<std::string> fam_ids_;
    std::vector<std::string> snp_ids_;
};

}  // namespace gelex
