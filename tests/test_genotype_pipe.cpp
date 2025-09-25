#include <catch2/catch_message.hpp>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <memory>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>

#include "gelex/data/genotype_pipe.h"
#include "gelex/data/sample_manager.h"
#include "gelex/error.h"

using Catch::Matchers::WithinAbs;

namespace test
{

inline std::filesystem::path replace_bed_extension(
    const std::string& bed_path,
    const std::string& new_extension)
{
    std::filesystem::path path(bed_path);
    path.replace_extension(new_extension);
    return path;
}

struct BinaryMatrixInfo
{
    int64_t num_samples;
    int64_t num_variants;
    std::vector<double> data;
};

struct SnpStatsInfo
{
    int64_t num_samples;
    int64_t num_variants;
    int64_t num_monomorphic;
    std::vector<int64_t> monomorphic_indices;
    std::vector<double> means;
    std::vector<double> stddevs;
};

BinaryMatrixInfo read_binary_matrix(
    const std::filesystem::path& file_path,
    int64_t expected_samples,
    int64_t expected_variants)
{
    std::ifstream file(file_path, std::ios::binary);
    REQUIRE(file.is_open());

    std::vector<double> data(expected_samples * expected_variants);
    const auto data_size
        = static_cast<std::streamsize>(data.size() * sizeof(double));
    file.read(reinterpret_cast<char*>(data.data()), data_size);

    REQUIRE(file.good());

    return BinaryMatrixInfo{
        .num_samples = expected_samples,
        .num_variants = expected_variants,
        .data = std::move(data)};
}

SnpStatsInfo read_snp_stats(const std::filesystem::path& file_path)
{
    std::ifstream file(file_path, std::ios::binary);
    REQUIRE(file.is_open());

    int64_t num_samples, num_variants, num_monomorphic;
    file.read(reinterpret_cast<char*>(&num_samples), sizeof(int64_t));
    file.read(reinterpret_cast<char*>(&num_variants), sizeof(int64_t));
    file.read(reinterpret_cast<char*>(&num_monomorphic), sizeof(int64_t));

    REQUIRE(file.good());

    std::vector<int64_t> monomorphic_indices(num_monomorphic);
    if (num_monomorphic > 0)
    {
        file.read(
            reinterpret_cast<char*>(monomorphic_indices.data()),
            static_cast<std::streamsize>(num_monomorphic * sizeof(int64_t)));
    }

    std::vector<double> means(num_variants);
    file.read(
        reinterpret_cast<char*>(means.data()),
        static_cast<std::streamsize>(num_variants * sizeof(double)));

    std::vector<double> stddevs(num_variants);
    file.read(
        reinterpret_cast<char*>(stddevs.data()),
        static_cast<std::streamsize>(num_variants * sizeof(double)));

    REQUIRE(file.good());

    return SnpStatsInfo{
        .num_samples = num_samples,
        .num_variants = num_variants,
        .num_monomorphic = num_monomorphic,
        .monomorphic_indices = std::move(monomorphic_indices),
        .means = std::move(means),
        .stddevs = std::move(stddevs)};
}
}  // namespace test

namespace test
{

class TestGenotypePipeManager
{
   public:
    explicit TestGenotypePipeManager(std::filesystem::path bed_file)
        : bed_file_{bed_file},
          fam_file_{bed_file.replace_extension(".fam")},
          bim_file_{bed_file.replace_extension(".bim")}
    {
        bed_ = *gelex::detail::open_file<std::ofstream>(
            bed_file_, std::ios::binary);
        bim_ = *gelex::detail::open_file<std::ofstream>(
            bim_file_, std::ios::out);
        fam_ = *gelex::detail::open_file<std::ofstream>(
            fam_file_, std::ios::out);
    }
    TestGenotypePipeManager(const TestGenotypePipeManager&) = delete;
    TestGenotypePipeManager(TestGenotypePipeManager&&) noexcept = default;
    TestGenotypePipeManager& operator=(const TestGenotypePipeManager&) = delete;
    TestGenotypePipeManager& operator=(TestGenotypePipeManager&&) noexcept
        = default;

    ~TestGenotypePipeManager()
    {
        std::filesystem::remove(bed_file_);
        std::filesystem::remove(fam_file_);
        std::filesystem::remove(bim_file_);

        // Clean up all possible output files
        std::filesystem::remove(bed_file_.replace_extension(".add.bmat"));
        std::filesystem::remove(bed_file_.replace_extension(".add.snpstats"));
        std::filesystem::remove(bed_file_.replace_extension(".dom.bmat"));
        std::filesystem::remove(bed_file_.replace_extension(".dom.snpstats"));

        // Clean up any files created by test sections that use different
        // prefixes
        std::filesystem::remove("test_genotype_pipe_mono.add.bmat");
        std::filesystem::remove("test_genotype_pipe_mono.add.snpstats");
        std::filesystem::remove("test_genotype_pipe_output.add.bmat");
        std::filesystem::remove("test_genotype_pipe_output.add.snpstats");
        std::filesystem::remove("test_genotype_pipe_single_sample.add.bmat");
        std::filesystem::remove(
            "test_genotype_pipe_single_sample.add.snpstats");
        std::filesystem::remove("test_genotype_pipe_single_variant.add.bmat");
        std::filesystem::remove(
            "test_genotype_pipe_single_variant.add.snpstats");
        std::filesystem::remove("test_genotype_pipe_file_exists.add.bmat");
        std::filesystem::remove("test_genotype_pipe_file_exists.add.snpstats");

        std::filesystem::remove("non_existent_file.bin");
        std::filesystem::remove("test_genotype_pipe_basic.add.snpstats");
        std::filesystem::remove("test_genotype_pipe_process.add.snpstats");
        std::filesystem::remove("test_output.bin");
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
    std::filesystem::path bed_file_;
    std::filesystem::path fam_file_;
    std::filesystem::path bim_file_;

    std::ofstream bed_;
    std::ofstream fam_;
    std::ofstream bim_;
};

}  // namespace test

TEST_CASE("GenotypePipe creation and basic functionality", "[genotype_pipe]")
{
    SECTION("Successful creation with valid files")
    {
        std::filesystem::path test_bed = "test_genotype_pipe_basic.bed";
        test::TestGenotypePipeManager pipe_manager(test_bed);

        std::vector<std::string> fids = {"fam1", "fam2", "fam3", "fam4"};
        std::vector<std::string> iids = {"ind1", "ind2", "ind3", "ind4"};
        pipe_manager.create(fids, iids, 5, 2);
        auto sample_manager
            = gelex::SampleManager::create(test_bed.replace_extension(".fam"))
                  .value();
        sample_manager.finalize();
        auto genotype_pipe = gelex::GenotypePipe::create(
            test_bed.replace_extension(".bed"),
            std::make_shared<gelex::SampleManager>(std::move(sample_manager)),
            test_bed.stem(),
            false);
        REQUIRE(genotype_pipe.has_value());
        REQUIRE(genotype_pipe->num_variants() == 5);

        auto result = genotype_pipe->process();
        REQUIRE(result.has_value());
        REQUIRE(genotype_pipe->num_samples() == 4);
        test::BinaryMatrixInfo bmat_info = test::read_binary_matrix(
            test::replace_bed_extension(test_bed, ".add.bmat"), 4, 5);

        Eigen::Map<const Eigen::MatrixXd> bmat(bmat_info.data.data(), 4, 5);

        Eigen::MatrixXd expected_bmat{
            {-1.22474487, 1.22474487, 0, 1.22474487, -1.22474487},
            {0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0.},
            {1.22474487, -1.22474487, 0, -1.22474487, 1.22474487}};
        REQUIRE(bmat.isApprox(expected_bmat, 1e-5));
    }

    SECTION("Creation with non-existent files")
    {
        auto genotype_pipe = gelex::GenotypePipe::create(
            "nonexistent", nullptr, "test_prefix", false);
        REQUIRE_FALSE(genotype_pipe.has_value());
        REQUIRE(genotype_pipe.error().code == gelex::ErrorCode::FileNotFound);
    }
}

TEST_CASE("GenotypePipe processing functionality", "[genotype_pipe]")
{
    std::filesystem::path test_bed = "test_genotype_pipe_process.bed";
    test::TestGenotypePipeManager pipe_manager(test_bed);

    std::vector<std::string> fids = {"fam1", "fam2", "fam3", "fam4"};
    std::vector<std::string> iids = {"ind1", "ind2", "ind3", "ind4"};
    pipe_manager.create(fids, iids, 6);

    auto sample_manager
        = gelex::SampleManager::create(test_bed.replace_extension(".fam"))
              .value();
    sample_manager.finalize();
    auto sample_ptr
        = std::make_shared<gelex::SampleManager>(std::move(sample_manager));

    auto genotype_pipe = gelex::GenotypePipe::create(
        test_bed.replace_extension(".bed"), sample_ptr, test_bed.stem(), false);

    auto result = genotype_pipe->process();
    REQUIRE(result.has_value());

    // Verify output files were created
    auto bmat_file = test::replace_bed_extension(test_bed, ".add.bmat");
    auto stats_file = test::replace_bed_extension(test_bed, ".add.snpstats");
    REQUIRE(std::filesystem::exists(bmat_file));
    REQUIRE(std::filesystem::exists(stats_file));

    // Verify file sizes are reasonable
    auto bmat_size = std::filesystem::file_size(bmat_file);
    auto stats_size = std::filesystem::file_size(stats_file);

    REQUIRE(bmat_size > 0);
    REQUIRE(stats_size > 0);
    test::BinaryMatrixInfo bmat_info
        = test::read_binary_matrix(bmat_file, 4, 6);
    test::SnpStatsInfo stats_info = test::read_snp_stats(stats_file);

    REQUIRE(bmat_info.num_samples == 4);
    REQUIRE(bmat_info.num_variants == 6);
    REQUIRE(bmat_info.data.size() == 24);

    REQUIRE(stats_info.num_samples == 4);
    REQUIRE(stats_info.num_variants == 6);
    REQUIRE(stats_info.num_monomorphic == 0);
    REQUIRE(stats_info.monomorphic_indices.empty());
    REQUIRE(stats_info.means.size() == 6);
    REQUIRE(stats_info.stddevs.size() == 6);

    std::filesystem::remove(bmat_file);
    std::filesystem::remove(stats_file);
    // Verify chunked processing produces same output
    // Create a new SampleManager for the second pipe
    auto genotype_pipe_chunk = gelex::GenotypePipe::create(
        test_bed, sample_ptr, test_bed.stem(), false);
    auto chunked_result = genotype_pipe_chunk->process(2);

    REQUIRE(std::filesystem::exists(bmat_file));
    REQUIRE(std::filesystem::exists(stats_file));

    test::BinaryMatrixInfo chunked_bmat_info
        = test::read_binary_matrix(bmat_file, 4, 6);
    test::SnpStatsInfo chunked_stats_info = test::read_snp_stats(stats_file);

    REQUIRE(chunked_bmat_info.num_samples == bmat_info.num_samples);
    REQUIRE(chunked_bmat_info.num_variants == bmat_info.num_variants);
    REQUIRE(chunked_bmat_info.data == bmat_info.data);
    REQUIRE(chunked_stats_info.num_samples == stats_info.num_samples);
    REQUIRE(chunked_stats_info.num_variants == stats_info.num_variants);
    REQUIRE(chunked_stats_info.num_monomorphic == stats_info.num_monomorphic);
    REQUIRE(
        chunked_stats_info.monomorphic_indices
        == stats_info.monomorphic_indices);
    REQUIRE(chunked_stats_info.means == stats_info.means);
    REQUIRE(chunked_stats_info.stddevs == stats_info.stddevs);

    // Verify file sizes match original
    auto chunked_bmat_size = std::filesystem::file_size(bmat_file);
    auto chunked_stats_size = std::filesystem::file_size(stats_file);

    REQUIRE(chunked_bmat_size == bmat_size);
    REQUIRE(chunked_stats_size == stats_size);
}

TEST_CASE("GenotypePipe monomorphic SNP detection", "[genotype_pipe]")
{
    std::filesystem::path test_bed = "test_genotype_pipe_mono.bed";
    test::TestGenotypePipeManager pipe_manager(test_bed);

    std::vector<std::string> fids = {"fam1", "fam2", "fam3", "fam4"};
    std::vector<std::string> iids = {"ind1", "ind2", "ind3", "ind4"};

    // Create data with one monomorphic SNP at index 2
    pipe_manager.create(fids, iids, 5, 2);

    auto sample_manager
        = gelex::SampleManager::create(test_bed.replace_extension(".fam"))
              .value();
    sample_manager.finalize();
    auto genotype_pipe = gelex::GenotypePipe::create(
        test_bed.replace_extension(".bed"),
        std::make_shared<gelex::SampleManager>(std::move(sample_manager)),
        test_bed.stem(),
        false);
    REQUIRE(genotype_pipe.has_value());

    auto result = genotype_pipe->process();
    REQUIRE(result.has_value());

    // Read the stats file to verify monomorphic detection
    std::ifstream stats_stream(
        "test_genotype_pipe_mono.add.snpstats", std::ios::binary);
    REQUIRE(stats_stream.is_open());

    int64_t num_samples = 0;
    int64_t num_variants = 0;
    int64_t num_monomorphic = 0;
    stats_stream.read(reinterpret_cast<char*>(&num_samples), sizeof(int64_t));
    stats_stream.read(reinterpret_cast<char*>(&num_variants), sizeof(int64_t));
    stats_stream.read(
        reinterpret_cast<char*>(&num_monomorphic), sizeof(int64_t));

    REQUIRE(num_samples == 4);
    REQUIRE(num_variants == 5);
    REQUIRE(num_monomorphic == 1);  // Should detect 1 monomorphic SNP

    // Read monomorphic indices
    std::vector<int64_t> monomorphic_indices(num_monomorphic);
    stats_stream.read(
        reinterpret_cast<char*>(monomorphic_indices.data()),
        static_cast<std::streamsize>(num_monomorphic * sizeof(int64_t)));

    REQUIRE(monomorphic_indices.size() == 1);
    REQUIRE(
        monomorphic_indices[0] == 2);  // SNP at index 2 should be monomorphic
}

TEST_CASE("GenotypePipe output verification", "[genotype_pipe]")
{
    std::filesystem::path test_bed = "test_genotype_pipe_output.bed";
    test::TestGenotypePipeManager pipe_manager(test_bed);

    std::vector<std::string> fids = {"fam1", "fam2"};
    std::vector<std::string> iids = {"ind1", "ind2"};
    pipe_manager.create(fids, iids, 3);

    auto sample_manager
        = gelex::SampleManager::create(test_bed.replace_extension(".fam"))
              .value();
    sample_manager.finalize();
    auto sample_ptr
        = std::make_shared<gelex::SampleManager>(std::move(sample_manager));
    auto genotype_pipe = gelex::GenotypePipe::create(
        test_bed.replace_extension(".bed"), sample_ptr, test_bed.stem(), false);

    REQUIRE(genotype_pipe.has_value());
    auto result = genotype_pipe->process();
    REQUIRE(result.has_value());

    // Load the original genotype data for comparison
    auto bed_sample_manager
        = gelex::SampleManager::create(test_bed.replace_extension(".fam"))
              .value();
    bed_sample_manager.finalize();
    auto bed_pipe = gelex::BedPipe::create(
        test_bed.replace_extension(".bed"), sample_ptr);
    REQUIRE(bed_pipe.has_value());

    auto original_matrix = bed_pipe->load();
    REQUIRE(original_matrix.has_value());

    auto bmat_path = test::replace_bed_extension(test_bed, ".add.bmat");
    std::ifstream bmat_stream(bmat_path, std::ios::binary);
    REQUIRE(bmat_stream.is_open());

    std::vector<double> processed_data(original_matrix->size());
    bmat_stream.read(
        reinterpret_cast<char*>(processed_data.data()),
        static_cast<std::streamsize>(original_matrix->size() * sizeof(double)));
    REQUIRE(processed_data.size() == original_matrix->size());

    // check snpstats
    auto stats_path = test::replace_bed_extension(test_bed, ".add.snpstats");
    std::ifstream stats_stream(stats_path, std::ios::binary);
    REQUIRE(stats_stream.is_open());

    int64_t num_samples = 0;
    int64_t num_variants = 0;
    int64_t num_monomorphic = 0;
    stats_stream.read(reinterpret_cast<char*>(&num_samples), sizeof(int64_t));
    stats_stream.read(reinterpret_cast<char*>(&num_variants), sizeof(int64_t));
    stats_stream.read(
        reinterpret_cast<char*>(&num_monomorphic), sizeof(int64_t));

    REQUIRE(num_samples == 2);
    REQUIRE(num_variants == 3);
    REQUIRE(num_monomorphic >= 0);
    REQUIRE(num_monomorphic <= 3);

    // Read the rest of the statistics
    if (num_monomorphic > 0)
    {
        std::vector<int64_t> monomorphic_indices(num_monomorphic);
        stats_stream.read(
            reinterpret_cast<char*>(monomorphic_indices.data()),
            static_cast<std::streamsize>(num_monomorphic * sizeof(int64_t)));

        for (auto idx : monomorphic_indices)
        {
            REQUIRE(idx >= 0);
            REQUIRE(idx < 3);
        }
    }

    std::vector<double> means(num_variants);
    stats_stream.read(
        reinterpret_cast<char*>(means.data()),
        static_cast<std::streamsize>(num_variants * sizeof(double)));

    std::vector<double> stddevs(num_variants);
    stats_stream.read(
        reinterpret_cast<char*>(stddevs.data()),
        static_cast<std::streamsize>(num_variants * sizeof(double)));

    // Verify statistics are reasonable
    for (int i = 0; i < num_variants; ++i)
    {
        REQUIRE(means[i] >= 0.0);
        REQUIRE(means[i] <= 2.0);  // Genotype values are 0,1,2
        REQUIRE(stddevs[i] >= 0.0);
    }
}

TEST_CASE("GenotypePipe edge cases", "[genotype_pipe]")
{
    SECTION("Process single sample")
    {
        std::filesystem::path test_bed = "test_genotype_pipe_single_sample.bed";
        test::TestGenotypePipeManager pipe_manager(test_bed);

        std::vector<std::string> fids = {"fam1"};
        std::vector<std::string> iids = {"ind1"};
        pipe_manager.create(fids, iids, 3);

        auto sample_manager
            = gelex::SampleManager::create(test_bed.replace_extension(".fam"))
                  .value();
        sample_manager.finalize();
        auto genotype_pipe = gelex::GenotypePipe::create(
            test_bed.replace_extension(".bed"),
            std::make_shared<gelex::SampleManager>(std::move(sample_manager)),
            test_bed.stem(),
            false);
        REQUIRE(genotype_pipe.has_value());

        auto result = genotype_pipe->process();
        REQUIRE(result.has_value());

        auto bmat_file = test::replace_bed_extension(test_bed, ".add.bmat");
        auto stats_file
            = test::replace_bed_extension(test_bed, ".add.snpstats");
        REQUIRE(std::filesystem::exists(bmat_file));
        REQUIRE(std::filesystem::exists(stats_file));
    }

    SECTION("Process single variant")
    {
        std::filesystem::path test_bed
            = "test_genotype_pipe_single_variant.bed";
        test::TestGenotypePipeManager pipe_manager(test_bed);

        std::vector<std::string> fids = {"fam1", "fam2"};
        std::vector<std::string> iids = {"ind1", "ind2"};
        pipe_manager.create(fids, iids, 1);

        auto sample_manager
            = gelex::SampleManager::create(test_bed.replace_extension(".fam"))
                  .value();
        sample_manager.finalize();
        auto genotype_pipe = gelex::GenotypePipe::create(
            test_bed.replace_extension(".bed"),
            std::make_shared<gelex::SampleManager>(std::move(sample_manager)),
            test_bed.stem(),
            false);
        REQUIRE(genotype_pipe.has_value());

        auto result = genotype_pipe->process();
        REQUIRE(result.has_value());

        auto bmat_file = test::replace_bed_extension(test_bed, ".add.bmat");
        auto stats_file
            = test::replace_bed_extension(test_bed, ".add.snpstats");
        REQUIRE(std::filesystem::exists(bmat_file));
        REQUIRE(std::filesystem::exists(stats_file));
    }
}

TEST_CASE("GenotypePipe file existence handling", "[genotype_pipe]")
{
    std::filesystem::path test_bed = "test_genotype_pipe_file_exists.bed";
    test::TestGenotypePipeManager pipe_manager(test_bed);

    std::vector<std::string> fids = {"fam1", "fam2"};
    std::vector<std::string> iids = {"ind1", "ind2"};
    pipe_manager.create(fids, iids, 3);

    auto sample_manager
        = gelex::SampleManager::create(test_bed.replace_extension(".fam"))
              .value();
    sample_manager.finalize();
    auto sample_ptr
        = std::make_shared<gelex::SampleManager>(std::move(sample_manager));

    SECTION("Skip processing when output files already exist")
    {
        // First run: create output files
        auto genotype_pipe1 = gelex::GenotypePipe::create(
            test_bed.replace_extension(".bed"),
            sample_ptr,
            test_bed.stem(),
            false);
        REQUIRE(genotype_pipe1.has_value());
        auto result1 = genotype_pipe1->process();
        REQUIRE(result1.has_value());

        auto bmat_file = test::replace_bed_extension(test_bed, ".add.bmat");
        auto stats_file
            = test::replace_bed_extension(test_bed, ".add.snpstats");
        REQUIRE(std::filesystem::exists(bmat_file));
        REQUIRE(std::filesystem::exists(stats_file));

        // Record file sizes and modification times
        auto original_bmat_size = std::filesystem::file_size(bmat_file);
        auto original_stats_size = std::filesystem::file_size(stats_file);
        auto original_bmat_time = std::filesystem::last_write_time(bmat_file);
        auto original_stats_time = std::filesystem::last_write_time(stats_file);

        // Second run: files exist, should skip processing
        auto genotype_pipe2 = gelex::GenotypePipe::create(
            test_bed.replace_extension(".bed"),
            sample_ptr,
            test_bed.stem(),
            false);

        REQUIRE_FALSE(genotype_pipe2.has_value());
        // Verify files were not modified (skipped processing)
        REQUIRE(std::filesystem::file_size(bmat_file) == original_bmat_size);
        REQUIRE(std::filesystem::file_size(stats_file) == original_stats_size);
        REQUIRE(
            std::filesystem::last_write_time(bmat_file) == original_bmat_time);
        REQUIRE(
            std::filesystem::last_write_time(stats_file)
            == original_stats_time);
    }

    SECTION("Process when output files don't exist")
    {
        auto bmat_file = test::replace_bed_extension(test_bed, ".add.bmat");
        auto stats_file
            = test::replace_bed_extension(test_bed, ".add.snpstats");

        // Ensure files don't exist initially
        if (std::filesystem::exists(bmat_file))
        {
            std::filesystem::remove(bmat_file);
        }
        if (std::filesystem::exists(stats_file))
        {
            std::filesystem::remove(stats_file);
        }

        REQUIRE_FALSE(std::filesystem::exists(bmat_file));
        REQUIRE_FALSE(std::filesystem::exists(stats_file));

        auto genotype_pipe = gelex::GenotypePipe::create(
            test_bed.replace_extension(".bed"),
            sample_ptr,
            test_bed.stem(),
            false);
        REQUIRE(genotype_pipe.has_value());
        auto result = genotype_pipe->process();
        REQUIRE(result.has_value());

        // Verify files were created
        REQUIRE(std::filesystem::exists(bmat_file));
        REQUIRE(std::filesystem::exists(stats_file));
        REQUIRE(std::filesystem::file_size(bmat_file) > 0);
        REQUIRE(std::filesystem::file_size(stats_file) > 0);
    }
}
