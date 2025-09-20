#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/data/bedpipe.h"
#include "gelex/error.h"

using Catch::Matchers::WithinAbs;

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
        bed_ = *gelex::detail::openfile<std::ofstream>(
            bed_file_, gelex::detail::file_type::binary);
        bim_ = *gelex::detail::openfile<std::ofstream>(bim_file_);
        fam_ = *gelex::detail::openfile<std::ofstream>(fam_file_);
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
            bim_ << std::format(
                "1\tsnp{}\t0\t{}\tG\tC\n", i + 1, ((i + 1) * 1000) + 1);
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

TEST_CASE("BedPipe creation and basic functionality", "[bed_pipe]")
{
    const std::string test_prefix = "test_bed_pipe_basic";
    test::TestBedManager bed_manager(test_prefix);

    SECTION("Successful creation with valid files")
    {
        std::vector<std::string> fids = {"fam1", "fam2", "fam3", "fam4"};
        std::vector<std::string> iids = {"ind1", "ind2", "ind3", "ind4"};
        bed_manager.create(fids, iids, 5);

        auto bed_pipe = gelex::BedPipe::create(test_prefix, false);
        REQUIRE(bed_pipe.has_value());

        REQUIRE(bed_pipe->sample_size() == 4);
        REQUIRE(bed_pipe->num_variants() == 5);

        const auto& sample_ids = bed_pipe->sample_map();
        REQUIRE(sample_ids.size() == 4);
        REQUIRE(sample_ids.contains("fam1_ind1"));
        REQUIRE(sample_ids.contains("fam2_ind2"));
        REQUIRE(sample_ids.contains("fam3_ind3"));
        REQUIRE(sample_ids.contains("fam4_ind4"));

        const auto& snp_ids = bed_pipe->snp_ids();
        REQUIRE(snp_ids.size() == 5);
        REQUIRE(snp_ids[0] == "snp1");
        REQUIRE(snp_ids[4] == "snp5");

        auto result = bed_pipe->load();
        REQUIRE(result.has_value());

        REQUIRE(result->rows() == 4);  // samples
        REQUIRE(result->cols() == 5);  // variants

        Eigen::MatrixXd expected{
            {0, 2, 0, 2, 0}, {1, 1, 1, 1, 1}, {1, 1, 1, 1, 1}, {2, 0, 2, 0, 2}};
        REQUIRE(result->isApprox(expected, 1e-10));
    }

    SECTION("Creation with IID-only mode")
    {
        std::vector<std::string> fids = {"fam1", "fam2"};
        std::vector<std::string> iids = {"ind1", "ind2"};
        bed_manager.create(fids, iids, 3);

        auto bed_pipe = gelex::BedPipe::create(test_prefix, true);
        REQUIRE(bed_pipe.has_value());
        REQUIRE(bed_pipe->sample_size() == 2);

        const auto& sample_ids = bed_pipe->sample_map();
        REQUIRE(sample_ids.size() == 2);
        REQUIRE(sample_ids.contains("ind1"));
        REQUIRE(sample_ids.contains("ind2"));
        REQUIRE(bed_pipe->snp_ids().size() == 3);
    }
}

TEST_CASE("BedPipe genotype access methods", "[bed_pipe]")
{
    const std::string test_prefix = "test_bed_pipe_genotypes";
    test::TestBedManager bed_manager(test_prefix);

    std::vector<std::string> fids = {"fam1", "fam2", "fam3", "fam4"};
    std::vector<std::string> iids = {"ind1", "ind2", "ind3", "ind4"};
    bed_manager.create(fids, iids, 5);

    auto bed_pipe = gelex::BedPipe::create(test_prefix);
    REQUIRE(bed_pipe.has_value());

    SECTION("get_genotypes for valid variant index")
    {
        auto genotypes = bed_pipe->get_genotypes(0);
        REQUIRE(genotypes.has_value());
        REQUIRE(genotypes->size() == 4);

        // Verify genotype values (based on our test pattern)
        // Pattern: 0b11100100 = 11 10 01 00 = 2.0, 1.0, NaN, 0.0

        REQUIRE(genotypes->isApprox(Eigen::VectorXd{{0, 1, 1, 2}}));
    }

    SECTION("get_genotype for valid indices")
    {
        auto genotype = bed_pipe->get_genotype(0, 0);
        REQUIRE(genotype.has_value());
        REQUIRE_THAT(*genotype, WithinAbs(0.0, 1e-10));

        genotype = bed_pipe->get_genotype(0, 1);
        REQUIRE(genotype.has_value());
        REQUIRE_THAT(*genotype, WithinAbs(1.0, 1e-10));

        genotype = bed_pipe->get_genotype(0, 2);
        REQUIRE(genotype.has_value());
        REQUIRE_THAT(*genotype, WithinAbs(1.0, 1e-10));

        genotype = bed_pipe->get_genotype(0, 3);
        REQUIRE(genotype.has_value());
        REQUIRE_THAT(*genotype, WithinAbs(2.0, 1e-10));
    }

    SECTION("get_sample_genotypes for valid sample index")
    {
        auto genotypes = bed_pipe->get_sample_genotypes(0);
        REQUIRE(genotypes.has_value());
        REQUIRE(genotypes->size() == 5);

        // Sample 0 should have consistent values across variants
        REQUIRE(genotypes->isApprox(Eigen::VectorXd{{0, 2, 0, 2, 0}}));
    }
}

TEST_CASE("BedPipe error handling", "[bed_pipe]")
{
    const std::string test_prefix = "test_bed_pipe_errors";
    test::TestBedManager bed_manager(test_prefix);

    std::vector<std::string> fids = {"fam1", "fam2"};
    std::vector<std::string> iids = {"ind1", "ind2"};
    bed_manager.create(fids, iids, 3);

    auto bed_pipe = gelex::BedPipe::create(test_prefix);
    REQUIRE(bed_pipe.has_value());

    SECTION("Invalid variant index")
    {
        auto genotypes = bed_pipe->get_genotypes(10);
        REQUIRE_FALSE(genotypes.has_value());
        REQUIRE(genotypes.error().code == gelex::ErrorCode::InvalidRange);

        auto genotype = bed_pipe->get_genotype(10, 0);
        REQUIRE_FALSE(genotype.has_value());
        REQUIRE(genotype.error().code == gelex::ErrorCode::InvalidRange);
    }

    SECTION("Invalid sample index")
    {
        auto genotype = bed_pipe->get_genotype(0, 10);
        REQUIRE_FALSE(genotype.has_value());
        REQUIRE(genotype.error().code == gelex::ErrorCode::InvalidRange);

        auto sample_genotypes = bed_pipe->get_sample_genotypes(10);
        REQUIRE_FALSE(sample_genotypes.has_value());
        REQUIRE(
            sample_genotypes.error().code == gelex::ErrorCode::InvalidRange);
    }

    SECTION("Invalid file prefix")
    {
        auto invalid_pipe = gelex::BedPipe::create("nonexistent_prefix");
        REQUIRE_FALSE(invalid_pipe.has_value());
        REQUIRE(invalid_pipe.error().code == gelex::ErrorCode::FileNotFound);
    }
}

TEST_CASE("BedPipe ID map validation and reordering", "[bed_pipe]")
{
    const std::string test_prefix = "test_bed_pipe_id_map";
    test::TestBedManager bed_manager(test_prefix);

    std::vector<std::string> fids = {"fam1", "fam2", "fam3", "fam4"};
    std::vector<std::string> iids = {"ind1", "ind2", "ind3", "ind4"};
    bed_manager.create(fids, iids, 5);

    auto bed_pipe = gelex::BedPipe::create(test_prefix, true);
    REQUIRE(bed_pipe.has_value());

    SECTION("Empty ID map validation")
    {
        std::unordered_map<std::string, Eigen::Index> empty_map;

        // Test direct validation
        auto validation = bed_pipe->load(empty_map);
        REQUIRE_FALSE(validation.has_value());
        REQUIRE(validation.error().code == gelex::ErrorCode::InvalidData);

        // Test chunk loading with empty map
        auto chunk_validation = bed_pipe->load_chunk(0, 2, empty_map);
        REQUIRE_FALSE(chunk_validation.has_value());
        REQUIRE(chunk_validation.error().code == gelex::ErrorCode::InvalidData);
    }

    SECTION("ID map with invalid sample IDs")
    {
        std::unordered_map<std::string, Eigen::Index> invalid_map
            = {{"ind1", 0}, {"nonexistent", 1}, {"ind3", 2}};

        auto validation = bed_pipe->load(invalid_map);
        REQUIRE_FALSE(validation.has_value());
        REQUIRE(validation.error().code == gelex::ErrorCode::InvalidData);
        REQUIRE(
            validation.error().message.find("nonexistent")
            != std::string::npos);
    }

    SECTION("Successful ID map reordering")
    {
        // Create a custom ID map with reordered samples
        std::unordered_map<std::string, Eigen::Index> custom_map = {
            {"ind4", 0},  // Original index 3 -> new index 0
            {"ind2", 1},  // Original index 1 -> new index 1
            {"ind1", 2}   // Original index 0 -> new index 2
        };

        auto reordered_matrix = bed_pipe->load(custom_map);
        REQUIRE(reordered_matrix.has_value());
        REQUIRE(reordered_matrix->rows() == 3);  // Only 3 samples in custom map
        REQUIRE(reordered_matrix->cols() == 5);  // All 5 variants

        // Load the original matrix for comparison
        auto original_matrix = bed_pipe->load();
        REQUIRE(original_matrix.has_value());
        REQUIRE(reordered_matrix->row(0).isApprox(original_matrix->row(3)));
        REQUIRE(reordered_matrix->row(1).isApprox(original_matrix->row(1)));
        REQUIRE(reordered_matrix->row(2).isApprox(original_matrix->row(0)));
    }

    SECTION("ID map reordering with chunk loading")
    {
        std::unordered_map<std::string, Eigen::Index> custom_map = {
            {"ind3", 0},  // Original index 2 -> new index 0
            {"ind1", 1}   // Original index 0 -> new index 1
        };

        auto reordered_chunk = bed_pipe->load_chunk(1, 4, custom_map);
        REQUIRE(reordered_chunk.has_value());
        REQUIRE(reordered_chunk->rows() == 2);  // 2 samples in custom map
        REQUIRE(reordered_chunk->cols() == 3);  // 3 variants (4-1)

        // Load original chunks for comparison
        auto original_chunk = bed_pipe->load_chunk(1, 4);
        REQUIRE(original_chunk.has_value());

        // Verify reordering is correct
        // ind3 (original row 2) should now be row 0
        REQUIRE(reordered_chunk->row(0).isApprox(original_chunk->row(2)));
        // ind1 (original row 0) should now be row 1
        REQUIRE(reordered_chunk->row(1).isApprox(original_chunk->row(0)));
    }

    SECTION("Optional ID map parameter with std::nullopt")
    {
        // Test that std::nullopt uses the original sample map
        auto matrix_with_nullopt = bed_pipe->load(std::nullopt);
        auto matrix_with_default = bed_pipe->load();

        REQUIRE(matrix_with_nullopt.has_value());
        REQUIRE(matrix_with_default.has_value());
        REQUIRE(matrix_with_nullopt->isApprox(*matrix_with_default));

        // Same test for load_chunk
        auto chunk_with_nullopt = bed_pipe->load_chunk(0, 2, std::nullopt);
        auto chunk_with_default = bed_pipe->load_chunk(0, 2);

        REQUIRE(chunk_with_nullopt.has_value());
        REQUIRE(chunk_with_default.has_value());
        REQUIRE(chunk_with_nullopt->isApprox(*chunk_with_default));
    }
}

TEST_CASE("BedPipe bulk loading", "[bed_pipe]")
{
    const std::string test_prefix = "test_bed_pipe_bulk";
    test::TestBedManager bed_manager(test_prefix);

    std::vector<std::string> fids = {"fam1", "fam2", "fam3", "fam4"};
    std::vector<std::string> iids = {"ind1", "ind2", "ind3", "ind4"};
    bed_manager.create(fids, iids, 10);

    auto bed_pipe = gelex::BedPipe::create(test_prefix);
    REQUIRE(bed_pipe.has_value());

    SECTION("Load entire genotype matrix")
    {
        auto matrix = bed_pipe->load();
        REQUIRE(matrix.has_value());

        REQUIRE(matrix->rows() == 4);   // samples
        REQUIRE(matrix->cols() == 10);  // variants

        auto first_variant = bed_pipe->get_genotypes(0);
        REQUIRE(first_variant.has_value());
        REQUIRE(matrix->col(0).isApprox(*first_variant));

        auto last_variant = bed_pipe->get_genotypes(9);
        REQUIRE(last_variant.has_value());
        REQUIRE(matrix->col(9).isApprox(*last_variant));
    }

    SECTION("Load chunk of genotype matrix")
    {
        auto chunk = bed_pipe->load_chunk(2, 6);
        REQUIRE(chunk.has_value());

        REQUIRE(chunk->rows() == 4);  // samples
        REQUIRE(chunk->cols() == 4);  // variants (6-2 = 4)

        for (Eigen::Index i = 0; i < 4; ++i)
        {
            auto variant = bed_pipe->get_genotypes(2 + i);
            REQUIRE(variant.has_value());
            REQUIRE(chunk->col(i).isApprox(*variant));
        }
    }

    SECTION("Invalid chunk range")
    {
        auto invalid_chunk = bed_pipe->load_chunk(5, 5);
        REQUIRE_FALSE(invalid_chunk.has_value());
        REQUIRE(invalid_chunk.error().code == gelex::ErrorCode::InvalidRange);

        invalid_chunk = bed_pipe->load_chunk(6, 5);
        REQUIRE_FALSE(invalid_chunk.has_value());
        REQUIRE(invalid_chunk.error().code == gelex::ErrorCode::InvalidRange);

        invalid_chunk = bed_pipe->load_chunk(8, 12);
        REQUIRE_FALSE(invalid_chunk.has_value());
        REQUIRE(invalid_chunk.error().code == gelex::ErrorCode::InvalidRange);
    }
}

TEST_CASE("BedPipe edge cases", "[bed_pipe]")
{
    const std::string test_prefix = "test_bed_pipe_edge";
    test::TestBedManager bed_manager(test_prefix);

    SECTION("Single sample")
    {
        std::vector<std::string> fids = {"fam1"};
        std::vector<std::string> iids = {"ind1"};
        bed_manager.create(fids, iids, 3);

        auto bed_pipe = gelex::BedPipe::create(test_prefix);
        REQUIRE(bed_pipe.has_value());
        REQUIRE(bed_pipe->sample_size() == 1);
        REQUIRE(bed_pipe->num_variants() == 3);

        auto genotypes = bed_pipe->get_genotypes(0);
        REQUIRE(genotypes.has_value());
        REQUIRE(genotypes->size() == 1);

        auto matrix = bed_pipe->load();
        REQUIRE(matrix.has_value());
        REQUIRE(matrix->rows() == 1);
        REQUIRE(matrix->cols() == 3);
    }

    SECTION("Single variant")
    {
        std::vector<std::string> fids = {"fam1", "fam2"};
        std::vector<std::string> iids = {"ind1", "ind2"};
        bed_manager.create(fids, iids, 1);

        auto bed_pipe = gelex::BedPipe::create(test_prefix);
        REQUIRE(bed_pipe.has_value());
        REQUIRE(bed_pipe->sample_size() == 2);
        REQUIRE(bed_pipe->num_variants() == 1);

        auto matrix = bed_pipe->load();
        REQUIRE(matrix.has_value());
        REQUIRE(matrix->rows() == 2);
        REQUIRE(matrix->cols() == 1);

        auto chunk = bed_pipe->load_chunk(0, 1);
        REQUIRE(chunk.has_value());
        REQUIRE(chunk->isApprox(*matrix));
    }

    SECTION("Monomorphic SNP")
    {
        std::vector<std::string> fids = {"fam1", "fam2", "fam3"};
        std::vector<std::string> iids = {"ind1", "ind2", "ind3"};
        bed_manager.create(fids, iids, 5, 2);

        auto bed_pipe = gelex::BedPipe::create(test_prefix);
        REQUIRE(bed_pipe.has_value());

        auto genotypes = bed_pipe->get_genotypes(2);
        REQUIRE(genotypes.has_value());

        for (Eigen::Index i = 0; i < genotypes->size(); ++i)
        {
            REQUIRE_THAT((*genotypes)(i), WithinAbs(0.0, 1e-10));
        }
    }
}
