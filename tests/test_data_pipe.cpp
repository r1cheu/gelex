#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <fstream>
#include <string>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/data/data_pipe.h"
#include "gelex/data/sample_manager.h"
#include "gelex/error.h"

using Catch::Matchers::WithinAbs;

class DataPipeTestFixture
{
   public:
    DataPipeTestFixture()
    {
        createValidFamFile();
        createValidPhenotypeFile();
        createValidQcovarFile();
        createValidCovarFile();
        createMalformedPhenotypeFile();
        createEmptyPhenotypeFile();
        createPhenotypeFileWithInvalidColumn();
    }

    DataPipeTestFixture(const DataPipeTestFixture&) = default;
    DataPipeTestFixture(DataPipeTestFixture&&) = delete;
    DataPipeTestFixture& operator=(const DataPipeTestFixture&) = default;
    DataPipeTestFixture& operator=(DataPipeTestFixture&&) = delete;
    ~DataPipeTestFixture()
    {
        std::remove("test_valid.fam");
        std::remove("test_valid.phe");
        std::remove("test_valid.qcovar");
        std::remove("test_valid.covar");
        std::remove("test_malformed.phe");
        std::remove("test_empty.phe");
        std::remove("test_invalid_column.phe");
    }

    static void createValidFamFile()
    {
        std::ofstream file("test_valid.fam");
        file << "FAM001 IND001 0 0 1 1\n"
             << "FAM001 IND002 0 0 2 1\n"
             << "FAM002 IND003 IND001 IND002 1 2\n"
             << "FAM003 IND004 0 0 1 -9\n"
             << "FAM004 IND005 IND003 IND004 2 1\n";
        file.flush();
    }

    static void createValidPhenotypeFile()
    {
        std::ofstream file("test_valid.phe");
        file << "FID\tIID\tsex\tseason\tday\tbwt\tloc\tdam\tT1\n"
             << "FAM001\tIND001\tMale\tWinter\t92\t1.2\tl32\tIND0921\t4.7658\n"
             << "FAM001\tIND002\tMale\tSpring\t88\t2.7\tl36\tIND0921\t12.4098\n"
             << "FAM002\tIND003\tMale\tSpring\t91\t1.0\tl17\tIND0968\t4.8545\n"
             << "FAM003\tIND004\tFemale\tAutumn\t82\t2.2\tl19\tIND1138\t36."
                "5418\n"
             << "FAM004\tIND005\tFemale\tSummer\t85\t1.8\tl25\tIND1201\t8."
                "9234\n";
    }

    static void createValidQcovarFile()
    {
        std::ofstream file("test_valid.qcovar");
        file << "FID\tIID\tage\tweight\theight\n"
             << "FAM001\tIND001\t25\t68.5\t175.2\n"
             << "FAM001\tIND002\t32\t72.1\t180.5\n"
             << "FAM002\tIND003\t28\t65.8\t172.8\n"
             << "FAM003\tIND004\t45\t78.3\t168.9\n"
             << "FAM004\tIND005\t36\t70.2\t177.1\n";
    }

    static void createValidCovarFile()
    {
        std::ofstream file("test_valid.covar");
        file << "FID\tIID\tsex\tlocation\tbatch\n"
             << "FAM001\tIND001\t1\t1\tA\n"
             << "FAM001\tIND002\t1\t2\tA\n"
             << "FAM002\tIND003\t1\t1\tB\n"
             << "FAM003\tIND004\t2\t3\tB\n"
             << "FAM004\tIND005\t2\t2\tA\n";
    }

    static void createMalformedPhenotypeFile()
    {
        std::ofstream file("test_malformed.phe");
        file
            << "FID\tIID\tsex\tseason\tday\tbwt\tloc\tdam\tT1\n"
            << "FAM001\tIND001\tMale\tWinter\t92\t1.2\tl32\tIND0921\t4.7658\n"
            << "FAM001\tIND002\tMale\tSpring\t88\t2.7\tl36\tIND0921\n";  // Missing
                                                                         // last
                                                                         // column
    }

    static void createEmptyPhenotypeFile()
    {
        std::ofstream file("test_empty.phe");
        // Empty file
    }

    static void createPhenotypeFileWithInvalidColumn()
    {
        std::ofstream file("test_invalid_column.phe");
        file << "FID\tIID\tsex\tseason\tday\tbwt\tloc\tdam\tT1\n"
             << "FAM001\tIND001\tMale\tWinter\t92\t1.2\tl32\tIND0921\t4.7658\n"
             << "FAM001\tIND002\tMale\tSpring\t88\t2.7\tl36\tIND0921\tinvalid_"
                "value\n";
    }
};

TEST_CASE_PERSISTENT_FIXTURE(
    DataPipeTestFixture,
    "DataPipe::create function",
    "[data_pipe][create]")
{
    SECTION("Valid configuration with all data types")
    {
        auto sample_manager
            = gelex::SampleManager::create("test_valid.fam", true);
        WARN(sample_manager.error().message);
        REQUIRE(sample_manager.has_value());

        gelex::DataPipe::Config config;
        config.phenotype_path = "test_valid.phe";
        config.phenotype_column = 8;
        config.qcovar_path = "test_valid.qcovar";
        config.covar_path = "test_valid.covar";
        config.iid_only = true;
        config.output_prefix = "test_output";

        auto result = gelex::DataPipe::create(
            config,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager.value())));
        REQUIRE(result.has_value());

        REQUIRE(result->has_phenotype());
        REQUIRE(result->has_fixed_effects());
        REQUIRE(result->phenotype().size() == 5);
        REQUIRE(result->fixed_effects().rows() == 5);
        // 3 qcovariates + (1 dummy vars for sex) + (2 dummy vars
        // for location) + (1 dummy var for batch)
        REQUIRE(result->fixed_effects().cols() == 1 + 3 + 1 + 2 + 1);
    }

    SECTION("Partial configuration - phenotype only")
    {
        auto sample_manager
            = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(sample_manager.has_value());

        gelex::DataPipe::Config config;
        config.phenotype_path = "test_valid.phe";
        config.phenotype_column = 8;
        config.iid_only = true;

        auto result = gelex::DataPipe::create(
            config,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager.value())));
        REQUIRE(result.has_value());

        REQUIRE(result->has_phenotype());
        REQUIRE(result->has_fixed_effects());
        REQUIRE(result->phenotype().size() == 5);
        REQUIRE(result->fixed_effects().cols() == 1);
    }

    SECTION("Partial configuration - qcovariates only")
    {
        auto sample_manager
            = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(sample_manager.has_value());

        gelex::DataPipe::Config config;
        config.phenotype_path = "test_valid.phe";
        config.phenotype_column = 8;
        config.qcovar_path = "test_valid.qcovar";
        config.iid_only = true;
        config.output_prefix = "test_output";

        auto result = gelex::DataPipe::create(
            config,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager.value())));
        REQUIRE(result.has_value());

        REQUIRE(result->has_phenotype());
        REQUIRE(result->has_fixed_effects());
        REQUIRE(result->phenotype().cols() == 1);
        REQUIRE(result->fixed_effects().rows() == 5);
        REQUIRE(result->fixed_effects().cols() == 4);  // 4 qcovariates
    }

    SECTION("Partial configuration - covariates only")
    {
        auto sample_manager
            = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(sample_manager.has_value());

        gelex::DataPipe::Config config;
        config.phenotype_path = "test_valid.phe";
        config.phenotype_column = 8;
        config.iid_only = true;
        config.output_prefix = "test_output";

        config.covar_path = "test_valid.covar";
        config.iid_only = true;

        auto result = gelex::DataPipe::create(
            config,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager.value())));
        REQUIRE(result.has_value());

        REQUIRE(result->has_phenotype());
        REQUIRE(result->has_fixed_effects());
        REQUIRE(result->phenotype().cols() == 1);
        REQUIRE(result->fixed_effects().rows() == 5);
        // intercept + (1 dummy vars for sex) + (2 dummy vars for location) + (1
        // dummy var for batch)
        REQUIRE(result->fixed_effects().cols() == 1 + 1 + 2 + 1);
    }

    SECTION("IID only vs full ID mode")
    {
        // Test IID only mode
        auto sample_manager_iid
            = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(sample_manager_iid.has_value());

        gelex::DataPipe::Config config_iid;
        config_iid.phenotype_path = "test_valid.phe";
        config_iid.phenotype_column = 8;
        config_iid.iid_only = true;

        auto result_iid = gelex::DataPipe::create(
            config_iid,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager_iid.value())));
        REQUIRE(result_iid.has_value());

        // Test full ID mode
        auto sample_manager_full
            = gelex::SampleManager::create("test_valid.fam", false);
        REQUIRE(sample_manager_full.has_value());

        gelex::DataPipe::Config config_full;
        config_full.phenotype_path = "test_valid.phe";
        config_full.phenotype_column = 8;
        config_full.iid_only = false;

        auto result_full = gelex::DataPipe::create(
            config_full,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager_full.value())));
        REQUIRE(result_full.has_value());

        // Both should have same number of samples
        REQUIRE(
            result_iid->phenotype().size() == result_full->phenotype().size());
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    DataPipeTestFixture,
    "DataPipe loading functionality",
    "[data_pipe][loading]")
{
    SECTION("Phenotype loading with different column indices")
    {
        auto sample_manager
            = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(sample_manager.has_value());

        // Test T1 column (index 8)
        gelex::DataPipe::Config config_t1;
        config_t1.phenotype_path = "test_valid.phe";
        config_t1.phenotype_column = 8;
        config_t1.iid_only = true;

        auto result_t1 = gelex::DataPipe::create(
            config_t1,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager.value())));
        REQUIRE(result_t1.has_value());
        REQUIRE(result_t1->phenotype_name() == "T1");

        // Test bwt column (index 5)
        auto sample_manager2
            = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(sample_manager2.has_value());

        gelex::DataPipe::Config config_bwt;
        config_bwt.phenotype_path = "test_valid.phe";
        config_bwt.phenotype_column = 5;
        config_bwt.iid_only = true;

        auto result_bwt = gelex::DataPipe::create(
            config_bwt,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager2.value())));
        REQUIRE(result_bwt.has_value());
        REQUIRE(result_bwt->phenotype_name() == "bwt");
    }

    SECTION("Qcovariates loading with multiple quantitative covariates")
    {
        auto sample_manager
            = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(sample_manager.has_value());

        gelex::DataPipe::Config config;
        config.phenotype_path = "test_valid.phe";
        config.phenotype_column = 8;
        config.qcovar_path = "test_valid.qcovar";
        config.iid_only = true;

        auto result = gelex::DataPipe::create(
            config,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager.value())));
        REQUIRE(result.has_value());

        REQUIRE(result->num_qcovariates() == 3);
        REQUIRE(result->qcovariate_names().size() == 3);
        REQUIRE(result->qcovariate_names()[0] == "age");
        REQUIRE(result->qcovariate_names()[1] == "weight");
        REQUIRE(result->qcovariate_names()[2] == "height");
    }

    SECTION("Covariates loading with categorical variables")
    {
        auto sample_manager
            = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(sample_manager.has_value());

        gelex::DataPipe::Config config;
        config.phenotype_path = "test_valid.phe";
        config.phenotype_column = 8;

        config.covar_path = "test_valid.covar";
        config.iid_only = true;

        auto result = gelex::DataPipe::create(
            config,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager.value())));
        REQUIRE(result.has_value());

        REQUIRE(result->num_covariates() == 3);
        REQUIRE(result->covariate_names().size() == 3);
        REQUIRE(result->covariate_names()[0] == "sex");
        REQUIRE(result->covariate_names()[1] == "location");
        REQUIRE(result->covariate_names()[2] == "batch");
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    DataPipeTestFixture,
    "DataPipe intersection and matrix conversion",
    "[data_pipe][intersection]")
{
    SECTION("Complete intersection with all samples")
    {
        auto sample_manager
            = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(sample_manager.has_value());

        gelex::DataPipe::Config config;
        config.phenotype_path = "test_valid.phe";
        config.phenotype_column = 8;
        config.qcovar_path = "test_valid.qcovar";
        config.covar_path = "test_valid.covar";
        config.iid_only = true;

        auto result = gelex::DataPipe::create(
            config,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager.value())));
        REQUIRE(result.has_value());

        REQUIRE(result->phenotype().size() == 5);
        REQUIRE(result->fixed_effects().rows() == 5);
        // intercept + 3 qcovariates + (1 dummy vars for sex) + (2 dummy vars
        // for location) + (1 dummy var for batch)
        REQUIRE(result->fixed_effects().cols() == 1 + 3 + 1 + 2 + 1);

        // Verify fixed effect names
        REQUIRE(result->fixed_effect_names().size() == 7);
        REQUIRE(result->fixed_effect_names()[1] == "age");
        REQUIRE(result->fixed_effect_names()[2] == "weight");
        REQUIRE(result->fixed_effect_names()[3] == "height");
        // The covariate names should include the original names, not the dummy
        // variable names
    }

    SECTION("Partial intersection with subset of samples")
    {
        // Create a FAM file with only 3 samples
        std::ofstream fam_file("test_partial.fam");
        fam_file << "FAM001 IND001 0 0 1 1\n"
                 << "FAM002 IND003 IND001 IND002 1 2\n"
                 << "FAM004 IND005 IND003 IND004 2 1\n";
        fam_file.close();

        auto sample_manager
            = gelex::SampleManager::create("test_partial.fam", true);
        REQUIRE(sample_manager.has_value());

        gelex::DataPipe::Config config;
        config.phenotype_path = "test_valid.phe";
        config.phenotype_column = 8;
        config.iid_only = true;

        auto result = gelex::DataPipe::create(
            config,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager.value())));
        REQUIRE(result.has_value());

        // Should only have 3 samples after intersection
        REQUIRE(result->phenotype().size() == 3);

        const auto& phenotype = result->phenotype();
        REQUIRE_THAT(phenotype(0), WithinAbs(4.7658, 1e-10));  // IND001
        REQUIRE_THAT(phenotype(1), WithinAbs(4.8545, 1e-10));  // IND003
        REQUIRE_THAT(phenotype(2), WithinAbs(8.9234, 1e-10));  // IND005

        std::remove("test_partial.fam");
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    DataPipeTestFixture,
    "DataPipe move semantics",
    "[data_pipe][move]")
{
    auto sample_manager = gelex::SampleManager::create("test_valid.fam", true);
    REQUIRE(sample_manager.has_value());

    gelex::DataPipe::Config config;
    config.phenotype_path = "test_valid.phe";
    config.phenotype_column = 8;
    config.qcovar_path = "test_valid.qcovar";
    config.iid_only = true;

    auto result = gelex::DataPipe::create(
        config,
        std::make_shared<gelex::SampleManager>(
            std::move(sample_manager.value())));
    REQUIRE(result.has_value());

    SECTION("take_phenotype method")
    {
        auto phenotype = std::move(*result).take_phenotype();
        REQUIRE(phenotype.size() == 5);
        REQUIRE_THAT(phenotype(0), WithinAbs(4.7658, 1e-10));

        // After move, original should be empty
        REQUIRE(result->phenotype().size() == 0);
    }

    SECTION("take_fixed_effects method")
    {
        auto fixed_effects = std::move(*result).take_fixed_effects();
        REQUIRE(fixed_effects.rows() == 5);
        REQUIRE(fixed_effects.cols() == 4);  // Intercept + 3 qcovariates

        // After move, original should be empty
        REQUIRE(result->fixed_effects().cols() == 0);
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    DataPipeTestFixture,
    "DataPipe error handling",
    "[data_pipe][error]")
{
    SECTION("File not found")
    {
        auto sample_manager
            = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(sample_manager.has_value());

        gelex::DataPipe::Config config;
        config.phenotype_path = "non_existent_file.phe";
        config.phenotype_column = 8;
        config.iid_only = true;

        auto result = gelex::DataPipe::create(
            config,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager.value())));
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::FileNotFound);
    }

    SECTION("Invalid column index")
    {
        auto sample_manager
            = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(sample_manager.has_value());

        gelex::DataPipe::Config config;
        config.phenotype_path = "test_valid.phe";
        config.phenotype_column = 1;  // Too low
        config.iid_only = true;

        auto result = gelex::DataPipe::create(
            config,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager.value())));
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InvalidRange);
    }

    SECTION("Malformed data file")
    {
        auto sample_manager
            = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(sample_manager.has_value());

        gelex::DataPipe::Config config;
        config.phenotype_path = "test_malformed.phe";
        config.phenotype_column = 8;
        config.iid_only = true;

        auto result = gelex::DataPipe::create(
            config,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager.value())));
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InconsistColumnCount);
    }

    SECTION("Empty phenotype file")
    {
        auto sample_manager
            = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(sample_manager.has_value());

        gelex::DataPipe::Config config;
        config.phenotype_path = "test_empty.phe";
        config.phenotype_column = 8;
        config.iid_only = true;

        auto result = gelex::DataPipe::create(
            config,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager.value())));
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InvalidFile);
    }

    SECTION("Invalid value in phenotype file")
    {
        auto sample_manager
            = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(sample_manager.has_value());

        gelex::DataPipe::Config config;
        config.phenotype_path = "test_invalid_column.phe";
        config.phenotype_column = 8;
        config.iid_only = true;

        auto result = gelex::DataPipe::create(
            config,
            std::make_shared<gelex::SampleManager>(
                std::move(sample_manager.value())));
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::NotNumber);
    }
}
