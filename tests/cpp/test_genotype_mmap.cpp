#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "../src/data/genotype_mmap.h"
#include "gelex/data/bed_io.h"
#include "gelex/data/io.h"

namespace test
{

class TestBedManager
{
   public:
    explicit TestBedManager(const std::string& prefix)
        : bed_file_{prefix + ".bed"},
          fam_file_{prefix + ".fam"},
          bim_file_{prefix + ".bim"}

    {
        bed_ = gelex::detail::open_or_throw<std::ofstream>(
            bed_file_, gelex::detail::file_type::binary);
        bim_ = gelex::detail::open_or_throw<std::ofstream>(bim_file_);
        fam_ = gelex::detail::open_or_throw<std::ofstream>(fam_file_);
    }
    TestBedManager(const TestBedManager&) = delete;
    TestBedManager(TestBedManager&&) noexcept = default;
    TestBedManager& operator=(const TestBedManager&) = delete;
    TestBedManager& operator=(TestBedManager&&) noexcept = default;

    ~TestBedManager()
    {
        std::filesystem::remove(bed_file_);
        std::filesystem::remove(fam_file_);
        std::filesystem::remove(bim_file_);
    }

    void create(
        std::span<const std::string> fids,
        std::span<const std::string> iids,
        int n_snps,
        std::optional<int> mono_snp_index = std::nullopt)
    {
        for (size_t i = 0; i < fids.size(); ++i)
        {
            fam_ << fids[i] << " " << iids[i] << " 0 0 0 -9\n";
        }
        fam_.flush();

        for (int i = 0; i < n_snps; ++i)
        {
            bim_ << "1 snp" << i + 1 << " 0 " << (i + 1) * 1000 << " A T\n";
        }
        bim_.flush();
        std::array<std::byte, 3> header
            = {std::byte{0x6C}, std::byte{0x1B}, std::byte{0x01}};
        bed_.write(reinterpret_cast<const char*>(header.data()), header.size());

        size_t bytes_per_snp = (fids.size() + 3) / 4;
        for (int snp = 0; snp < n_snps; ++snp)
        {
            for (size_t byte_idx = 0; byte_idx < bytes_per_snp; ++byte_idx)
            {
                if (mono_snp_index && snp == mono_snp_index.value())
                {
                    // Create monomorphic SNP (all 00 genotypes)
                    unsigned char genotype_byte = 0b00000000;
                    bed_.write(
                        reinterpret_cast<const char*>(&genotype_byte), 1);
                }
                else
                {
                    // Create polymorphic SNP with variation
                    unsigned char genotype_byte
                        = 0b11100100;  // 11 10 01 00 in binary
                    if (snp % 2 == 1)
                    {
                        genotype_byte = 0b00011011;  // 00 01 10 11 in binary
                    }
                    bed_.write(
                        reinterpret_cast<const char*>(&genotype_byte), 1);
                }
            }
        }
        bed_.flush();
    }

   private:
    std::string bed_file_;
    std::string fam_file_;
    std::string bim_file_;

    std::ofstream bed_;
    std::ofstream fam_;
    std::ofstream bim_;
};
}  // namespace test

TEST_CASE("create_genotype_binary basic functionality", "[genotype_mmap]")
{
    const std::string test_prefix = "test_genotype_binary";
    test::TestBedManager bed_manager(test_prefix);
    const std::string bin_file = test_prefix + ".add.bin";
    const std::string meta_file = test_prefix + ".add.meta";

    std::filesystem::remove(bin_file);
    std::filesystem::remove(meta_file);

    SECTION("Valid input files")
    {
        // Create test data files
        std::vector<std::string> fids = {"fam1", "fam2", "fam3", "fam4"};
        std::vector<std::string> iids = {"ind1", "ind2", "ind3", "ind4"};

        bed_manager.create(fids, iids, 5);

        std::vector<std::string> p_list = {"ind1", "ind2", "ind3"};
        std::vector<std::string> g_list = {"ind1", "ind2", "ind3", "ind4"};

        REQUIRE_NOTHROW(
            gelex::detail::create_genotype_binary(
                test_prefix, false, p_list, g_list, true));

        // Check that output files were created
        REQUIRE(std::filesystem::exists(bin_file));
        REQUIRE(std::filesystem::exists(meta_file));
    }

    SECTION("Invalid bed file throws exception")
    {
        std::vector<std::string> p_list = {"ind1"};
        std::vector<std::string> g_list = {"ind1"};

        REQUIRE_THROWS_AS(
            gelex::detail::create_genotype_binary(
                "nonexistent", false, p_list, g_list, true),
            std::runtime_error);
    }

    SECTION("Phenotype ID not in genotype throws exception")
    {
        std::vector<std::string> fids = {"fam1", "fam2"};
        std::vector<std::string> iids = {"ind1", "ind2"};
        bed_manager.create(fids, iids, 2);

        std::vector<std::string> p_list
            = {"ind1", "ind_missing"};  // ind_missing not in fam
        std::vector<std::string> g_list = {"ind1", "ind2"};

        REQUIRE_THROWS_AS(
            gelex::detail::create_genotype_binary(
                test_prefix, false, p_list, g_list, true),
            std::runtime_error);
    }

    std::filesystem::remove(bin_file);
    std::filesystem::remove(meta_file);
}

TEST_CASE("create_genotype_binary with family IDs", "[genotype_mmap]")
{
    const std::string test_prefix = "test_genotype_binary_fam";
    test::TestBedManager bed_manager(test_prefix);
    const std::string bin_file = test_prefix + ".add.bin";
    const std::string meta_file = test_prefix + ".add.meta";

    // Cleanup existing files
    std::filesystem::remove(bin_file);
    std::filesystem::remove(meta_file);

    SECTION("Using FID_IID format")
    {
        std::vector<std::string> fids = {"fam1", "fam2", "fam3"};
        std::vector<std::string> iids = {"ind1", "ind2", "ind3"};
        bed_manager.create(fids, iids, 3);

        std::vector<std::string> p_list = {"fam1_ind1", "fam2_ind2"};
        std::vector<std::string> g_list
            = {"fam1_ind1", "fam2_ind2", "fam3_ind3"};

        REQUIRE_NOTHROW(
            gelex::detail::create_genotype_binary(
                test_prefix, false, p_list, g_list, false));

        REQUIRE(std::filesystem::exists(bin_file));
        REQUIRE(std::filesystem::exists(meta_file));
    }

    // Cleanup
    std::filesystem::remove(bin_file);
    std::filesystem::remove(meta_file);
}

TEST_CASE("GenotypeMap initialization", "[genotype_mmap]")
{
    const std::string test_prefix = "test_genotype_map";
    const std::string bin_file = test_prefix + ".add.bin";
    const std::string meta_file = test_prefix + ".add.meta";

    test::TestBedManager bed_manager(test_prefix);

    // Cleanup existing files
    std::filesystem::remove(bin_file);
    std::filesystem::remove(meta_file);

    SECTION("Valid binary file")
    {
        // First create a binary file using create_genotype_binary
        std::vector<std::string> fids = {"fam1", "fam2", "fam3", "fam4"};
        std::vector<std::string> iids = {"ind1", "ind2", "ind3", "ind4"};
        bed_manager.create(fids, iids, 5);

        std::vector<std::string> p_list = {"ind1", "ind2", "ind3"};
        std::vector<std::string> g_list = {"ind1", "ind2", "ind3", "ind4"};

        gelex::detail::create_genotype_binary(
            test_prefix, false, p_list, g_list, true);

        // Now test GenotypeMap
        gelex::detail::GenotypeMap gmap(bin_file);

        REQUIRE(gmap.mat.rows() == 3);  // p_list size
        REQUIRE(gmap.mat.cols() == 5);  // number of SNPs

        REQUIRE_NOTHROW(gmap.mat(0, 0));
        REQUIRE_NOTHROW(gmap.mat(2, 4));
    }

    SECTION("Invalid binary file throws exception")
    {
        REQUIRE_THROWS_AS(
            gelex::detail::GenotypeMap("nonexistent.bin"), std::runtime_error);
    }

    SECTION("File too small throws exception")
    {
        std::ofstream small_file(bin_file, std::ios::binary);
        small_file.write("small", 5);
        small_file.close();

        REQUIRE_THROWS_AS(
            gelex::detail::GenotypeMap(bin_file), std::runtime_error);
    }

    SECTION("Invalid header dimensions throws exception")
    {
        // Create a file with invalid dimensions in header
        std::ofstream bad_file(bin_file, std::ios::binary);
        int64_t bad_rows = -1;
        int64_t bad_cols = 5;
        bad_file.write(
            reinterpret_cast<const char*>(&bad_rows), sizeof(int64_t));
        bad_file.write(
            reinterpret_cast<const char*>(&bad_cols), sizeof(int64_t));
        bad_file.close();

        REQUIRE_THROWS_AS(
            gelex::detail::GenotypeMap(bin_file), std::runtime_error);
    }

    // Cleanup
    std::filesystem::remove(bin_file);
    std::filesystem::remove(meta_file);
}

TEST_CASE("GenotypeMap matrix data integrity", "[genotype_mmap]")
{
    const std::string test_prefix = "test_genotype_map_data";
    const std::string bin_file = test_prefix + ".add.bin";
    const std::string meta_file = test_prefix + ".add.meta";

    test::TestBedManager bed_manager(test_prefix);
    // Cleanup existing files
    std::filesystem::remove(bin_file);
    std::filesystem::remove(meta_file);

    SECTION("Matrix data consistency")
    {
        std::vector<std::string> fids = {"fam1", "fam2"};
        std::vector<std::string> iids = {"ind1", "ind2"};
        bed_manager.create(fids, iids, 3);

        std::vector<std::string> p_list = {"ind1", "ind2"};
        std::vector<std::string> g_list = {"ind1", "ind2"};

        gelex::detail::create_genotype_binary(
            test_prefix, false, p_list, g_list, true);

        gelex::detail::GenotypeMap gmap(bin_file);

        // Test matrix dimensions
        REQUIRE(gmap.mat.rows() == 2);
        REQUIRE(gmap.mat.cols() == 3);

        Eigen::MatrixXd expected(2, 3);
        expected << 0.7071, -0.7071, 0.7071, -0.7071, 0.7071, -0.7071;

        // Check matrix values with tolerance
        REQUIRE(expected.isApprox(gmap.mat, 1e-4));
    }

    SECTION("Matrix data consistency, change g order")
    {
        std::vector<std::string> fids = {"fam1", "fam2"};
        std::vector<std::string> iids = {"ind1", "ind2"};
        bed_manager.create(fids, iids, 3);

        std::vector<std::string> p_list = {"ind1", "ind2"};
        std::vector<std::string> g_list = {"ind2", "ind1"};

        gelex::detail::create_genotype_binary(
            test_prefix, false, p_list, g_list, true);

        gelex::detail::GenotypeMap gmap(bin_file);

        // Test matrix dimensions
        REQUIRE(gmap.mat.rows() == 2);
        REQUIRE(gmap.mat.cols() == 3);

        Eigen::MatrixXd expected{
            {0.7071, -0.7071, 0.7071}, {-0.7071, 0.7071, -0.7071}};

        // Check matrix values with tolerance
        REQUIRE(expected.isApprox(gmap.mat, 1e-4));
    }

    SECTION("Matrix data consistency, change g order")
    {
        std::vector<std::string> fids = {"fam1", "fam2"};
        std::vector<std::string> iids = {"ind1", "ind2"};
        bed_manager.create(fids, iids, 3);

        std::vector<std::string> p_list = {"ind2", "ind1"};
        std::vector<std::string> g_list = {"ind2", "ind1"};

        gelex::detail::create_genotype_binary(
            test_prefix, false, p_list, g_list, true);

        gelex::detail::GenotypeMap gmap(bin_file);

        // Test matrix dimensions
        REQUIRE(gmap.mat.rows() == 2);
        REQUIRE(gmap.mat.cols() == 3);

        Eigen::MatrixXd expected{
            {-0.7071, 0.7071, -0.7071}, {0.7071, -0.7071, 0.7071}};

        // Check matrix values with tolerance
        REQUIRE(expected.isApprox(gmap.mat, 1e-4));
    }

    // Cleanup
    std::filesystem::remove(bin_file);
    std::filesystem::remove(meta_file);
}

TEST_CASE("create_genotype_binary with monomorphic SNPs", "[genotype_mmap]")
{
    const std::string test_prefix = "test_genotype_binary_mono";
    const std::string bin_file = test_prefix + ".add.bin";
    const std::string meta_file = test_prefix + ".add.meta";

    test::TestBedManager bed_manager(test_prefix);

    // Cleanup existing files
    std::filesystem::remove(bin_file);
    std::filesystem::remove(meta_file);

    SECTION("Detect monomorphic SNPs")
    {
        // Create test data with one monomorphic SNP
        std::vector<std::string> fids = {"fam1", "fam2", "fam3", "fam4"};
        std::vector<std::string> iids = {"ind1", "ind2", "ind3", "ind4"};
        bed_manager.create(fids, iids, 5, 2);

        std::vector<std::string> p_list = {"ind1", "ind2", "ind3"};
        std::vector<std::string> g_list = {"ind1", "ind2", "ind3", "ind4"};

        REQUIRE_NOTHROW(
            gelex::detail::create_genotype_binary(
                test_prefix, false, p_list, g_list, true));

        // Check that output files were created
        REQUIRE(std::filesystem::exists(bin_file));
        REQUIRE(std::filesystem::exists(meta_file));

        // Test that we can read the mono indices back
        gelex::detail::GenotypeMap gmap(bin_file);
        REQUIRE(gmap.mat.rows() == 3);
        REQUIRE(gmap.mat.cols() == 5);

        // The monomorphic SNP should have zero variance
        Eigen::VectorXd snp2 = gmap.mat.col(2);
        double variance
            = (snp2.array() - snp2.mean()).square().sum() / (snp2.size() - 1);
        REQUIRE(variance < 1e-9);  // Should be effectively zero
    }

    // Cleanup
    std::filesystem::remove(bin_file);
    std::filesystem::remove(meta_file);
}

TEST_CASE("create_genotype_binary with dominance", "[genotype_mmap]")
{
    const std::string test_prefix = "test_genotype_binary_dom";
    const std::string add_bin = test_prefix + ".add.bin";
    const std::string add_meta = test_prefix + ".add.meta";
    const std::string dom_bin = test_prefix + ".dom.bin";
    const std::string dom_meta = test_prefix + ".dom.meta";

    test::TestBedManager bed_manager(test_prefix);
    // Cleanup existing files
    std::filesystem::remove(add_bin);
    std::filesystem::remove(add_meta);
    std::filesystem::remove(dom_bin);
    std::filesystem::remove(dom_meta);

    SECTION("Generate both additive and dominance files")
    {
        // Create test data
        std::vector<std::string> fids = {"fam1", "fam2", "fam3"};
        std::vector<std::string> iids = {"ind1", "ind2", "ind3"};
        bed_manager.create(fids, iids, 3);

        std::vector<std::string> p_list = {"ind1", "ind2"};
        std::vector<std::string> g_list = {"ind1", "ind2", "ind3"};

        REQUIRE_NOTHROW(
            gelex::detail::create_genotype_binary(
                test_prefix, true, p_list, g_list, true));

        // Check that both additive and dominance files were created
        REQUIRE(std::filesystem::exists(add_bin));
        REQUIRE(std::filesystem::exists(add_meta));
        REQUIRE(std::filesystem::exists(dom_bin));
        REQUIRE(std::filesystem::exists(dom_meta));

        // Test that both genotype maps can be loaded
        gelex::detail::GenotypeMap add_map(add_bin);
        gelex::detail::GenotypeMap dom_map(dom_bin);

        REQUIRE(add_map.mat.rows() == 2);
        REQUIRE(add_map.mat.cols() == 3);
        REQUIRE(dom_map.mat.rows() == 2);
        REQUIRE(dom_map.mat.cols() == 3);

        // Dominance should have different values than additive
        REQUIRE_FALSE(add_map.mat.isApprox(dom_map.mat));
    }

    // Cleanup
    std::filesystem::remove(add_bin);
    std::filesystem::remove(add_meta);
    std::filesystem::remove(dom_bin);
    std::filesystem::remove(dom_meta);
}
