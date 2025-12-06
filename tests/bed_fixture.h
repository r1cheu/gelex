#ifndef GELEX_TEST_BED_FIXTURE_H
#define GELEX_TEST_BED_FIXTURE_H

#include <cstddef>
#include <filesystem>
#include <random>
#include <span>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <utility>

#include "file_fixture.h"

namespace gelex::test
{

class BedFixture
{
   public:
    BedFixture();
    ~BedFixture() = default;

    BedFixture(const BedFixture&) = delete;
    BedFixture& operator=(const BedFixture&) = delete;
    BedFixture(BedFixture&&) noexcept = default;
    BedFixture& operator=(BedFixture&&) noexcept = default;

    std::pair<std::filesystem::path, Eigen::MatrixXd> create_bed_files(
        Eigen::Index num_samples,
        Eigen::Index num_snps,
        double missing_rate = 0.0,
        double maf_min = 0.05,
        double maf_max = 0.5,
        uint64_t seed = std::random_device{}());

    [[nodiscard]] std::filesystem::path get_prefix() const noexcept
    {
        return current_prefix_;
    }

    [[nodiscard]] FileFixture& get_file_fixture() noexcept
    {
        return file_fixture_;
    }

   private:
    FileFixture file_fixture_;
    std::filesystem::path current_prefix_;
    std::mt19937_64 rng_;

    static std::vector<std::byte> encode_variant(
        const Eigen::VectorXd& variant);

    void write_bed_file(const Eigen::MatrixXd& genotypes);

    void write_bim_file(
        Eigen::Index num_snps,
        std::span<const std::string> snp_ids,
        std::span<const std::string> chromosomes);

    void write_fam_file(
        Eigen::Index num_samples,
        std::span<const std::string> sample_ids);

    static std::pair<char, char> generate_random_alleles(std::mt19937_64& rng);

    static std::vector<std::string> generate_random_sample_ids(
        Eigen::Index num_samples,
        std::mt19937_64& rng);

    static std::vector<std::string> generate_random_snp_ids(
        Eigen::Index num_snps,
        std::mt19937_64& rng);

    static std::vector<std::string> generate_random_chromosomes(
        Eigen::Index num_snps,
        std::mt19937_64& rng);
};

}  // namespace gelex::test

#endif  // GELEX_TEST_BED_FIXTURE_H
