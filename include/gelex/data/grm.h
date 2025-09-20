#pragma once
#include <memory>
#include <string>
#include <unordered_set>

#include <Eigen/Core>

#include "gelex/data/bedpipe.h"

namespace gelex
{

class GRM
{
   public:
    explicit GRM(
        std::string_view bed_file,
        Eigen::Index chunk_size = 10000,
        const std::vector<std::string>& target_order = {});
    GRM(const GRM&) = delete;
    GRM(GRM&&) noexcept = default;
    GRM& operator=(const GRM&) = delete;
    GRM& operator=(GRM&&) noexcept = default;

    ~GRM() = default;

    Eigen::MatrixXd compute(bool add = true);
    Eigen::VectorXd p_major() const noexcept { return p_major_; }
    double scale_factor() const { return scale_factor_; }

   private:
    double scale_factor_{1.0};
    Eigen::Index chunk_size_{10000};

    Eigen::VectorXd p_major_;
    std::unique_ptr<BedPipe> bed_;
};

class CrossGRM
{
   public:
    explicit CrossGRM(
        std::string_view train_bed,
        Eigen::VectorXd p_major,
        double scale_factor,
        Eigen::Index chunk_size = 10000,
        const std::vector<std::string>& target_order = {});

    CrossGRM(const CrossGRM&) = delete;
    CrossGRM(CrossGRM&&) noexcept = default;
    CrossGRM& operator=(const CrossGRM&) = delete;
    CrossGRM& operator=(CrossGRM&&) noexcept = default;

    ~CrossGRM() = default;

    Eigen::MatrixXd compute(std::string_view test_bed, bool add = true);

   private:
    void check_snp_consistency(const BedPipe& test_bed) const;

    std::unique_ptr<BedPipe> bed_;
    Eigen::VectorXd p_major_;
    double scale_factor_;

    Eigen::Index chunk_size_;
};

/**
 * @brief Compute the allele frequencies for each SNP in the genotype matrix.
 *
 * This function calculates the mean of each column (SNP) in the genotype matrix
 * and divides by 2 to obtain the allele frequency. The result is a vector where
 * each element represents the frequency of the reference allele for the
 * corresponding SNP.
 *
 * @param genotype The genotype matrix where rows represent individuals and
 * columns represent SNPs
 * @return arma::dvec Vector of allele frequencies (mean/2) for each SNP
 */
Eigen::VectorXd compute_p_major(
    const Eigen::Ref<const Eigen::MatrixXd>& genotype);

void code_method_varden(
    const Eigen::Ref<const Eigen::VectorXd>& p_major,
    Eigen::Ref<Eigen::MatrixXd> genotype,
    bool add = true);

}  // namespace gelex
