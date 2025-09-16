#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "gelex/data/io.h"

// Helper class to manage temporary files for tests.
// Creates a directory and files in the constructor, and cleans them up in the
// destructor.
class DataReaderTestFixture
{
   public:
    DataReaderTestFixture()
    {
        // Create a unique directory for test files to avoid conflicts.
        test_dir = "datareader_test_files";
        std::filesystem::create_directory(test_dir);

        // Define file paths
        pheno_path_ = test_dir / "pheno.txt";
        fam_path_ = test_dir / "fam.txt";
        qcovar_path_ = test_dir / "qcovar.txt";
        covar_path_ = test_dir / "covar.txt";
    }

    ~DataReaderTestFixture()
    {
        // Cleanup: remove the directory and all its contents.
        std::filesystem::remove_all(test_dir);
    }

    // Method to create files for the corrected "happy path" scenario.
    void create_happy_path_files()
    {
        // .fam file with 4 individuals
        std::ofstream fam_file(fam_path_);
        fam_file << "FID1 IID1 0 0 1 1\n"
                 << "FID2 IID2 0 0 1 1\n"
                 << "FID3 IID3 0 0 1 1\n"
                 << "FID4 IID4 0 0 1 1\n";
        fam_file.close();

        // phenotype file, IID4 is not present.
        std::ofstream pheno_file(pheno_path_);
        pheno_file << "FID\tIID\tPHENO1\n"
                   << "FID1\tIID1\t10.5\n"
                   << "FID2\tIID2\t-20.2\n"
                   << "FID3\tIID3\t30.0\n"
                   << "FID4\tIID4\tNA\n";
        pheno_file.close();

        // quantitative covariate file, IID4 is not present.
        std::ofstream qcovar_file(qcovar_path_);
        qcovar_file << "FID\tIID\tQC1\tQC2\n"
                    << "FID1\tIID1\t1.1\t2.2\n"
                    << "FID2\tIID2\t3.3\t4.4\n"
                    << "FID3\tIID3\t5.5\t6.6\n"
                    << "FID4\tIID4\t7.7\t8.8\n";
        qcovar_file.close();

        // categorical covariate file.
        // For the intersected individuals (IID1, IID2, IID3):
        // C1 has two levels: 'A', 'B'.
        // C2 has three levels: 'X', 'Y', 'Z'.
        std::ofstream covar_file(covar_path_);
        covar_file << "FID\tIID\tC1\tC2\n"
                   << "FID1\tIID1\tA\tX\n"
                   << "FID2\tIID2\tB\tY\n"
                   << "FID3\tIID3\tA\tZ\n"
                   << "FID4\tIID4\tC\tX\n"   // IID4 will be excluded
                   << "FID6\tIID6\tC\tX\n";  // IID6 will be excluded
        covar_file.close();
    }

    // Method to create files for an invalid header scenario.
    void create_invalid_header_files()
    {
        // Valid .fam file
        std::ofstream fam_file(fam_path_);
        fam_file << "FID1 IID1 0 0 1 1\n";
        fam_file.close();

        // Phenotype file with an invalid header
        std::ofstream pheno_file(pheno_path_);
        pheno_file << "BAD_FID\tBAD_IID\tPHENO1\n"  // Invalid header
                   << "FID1\tIID1\t10.5\n";
        pheno_file.close();
    }

   protected:
    std::filesystem::path test_dir;
    std::filesystem::path pheno_path_;
    std::filesystem::path fam_path_;
    std::filesystem::path qcovar_path_;
    std::filesystem::path covar_path_;
};

TEST_CASE_METHOD(DataReaderTestFixture, "DataReader Class Unit Tests")
{
    SECTION("Happy Path: Covariates with 2 and 3 levels respectively")
    {
        // Description:
        // 此测试确保当所有输入文件有效，且交集后的样本中，一个分类协变量有2个水平，
        // 另一个有3个水平时，程序能正确进行独热编码并构建固定效应矩阵。
        create_happy_path_files();

        // 根据文件设置，IID1, IID2, IID3 在所有四个文件中都存在且数据有效。
        auto reader = gelex::DataReader::Create(
            pheno_path_, fam_path_, qcovar_path_, covar_path_, 3, true);

        // 交集后，预期 final_ids 按字典序排序为 {"IID1", "IID2", "IID3"}。
        const auto& final_ids = reader.final_ids();
        REQUIRE(final_ids.size() == 3);
        REQUIRE(final_ids[0] == "IID1");
        REQUIRE(final_ids[1] == "IID2");
        REQUIRE(final_ids[2] == "IID3");

        // 验证表型向量
        const auto& phenotype = reader.phenotype();
        REQUIRE(phenotype.size() == 3);
        CHECK(phenotype(0) == 10.5);   // IID1
        CHECK(phenotype(1) == -20.2);  // IID2
        CHECK(phenotype(2) == 30.0);   // IID3

        // 验证固定效应矩阵
        const auto& fixed = reader.fixed();
        // 行数 = 3 (个体数)
        // 列数 = 1 (截距) + 2 (qcovar: QC1, QC2) + 1 (C1有2个水平) + 2
        // (C2有3个水平) = 6
        //
        std::cout << fixed;
        REQUIRE(fixed.rows() == 3);
        REQUIRE(fixed.cols() == 6);

        // 预期矩阵内容: [截距, QC1, QC2, C1_B, C2_Y, C2_Z]
        // 编码规则:
        // C1 水平排序为 A, B。A 作为参考(0)，B 编码为 1。 -> 1列
        // C2 水平排序为 X, Y, Z。X 作为参考(0,0)，Y 编码为(1,0)，Z
        // 编码为(0,1)。 -> 2列
        Eigen::MatrixXd expected_fixed(3, 6);
        //                intercept, QC1, QC2, C1_B, C2_Y, C2_Z
        expected_fixed << 1.0, 1.1, 2.2, 0.0, 0.0, 0.0,  // IID1 (C1=A, C2=X)
            1.0, 3.3, 4.4, 1.0, 1.0, 0.0,                // IID2 (C1=B, C2=Y)
            1.0, 5.5, 6.6, 0.0, 0.0, 1.0;                // IID3 (C1=A, C2=Z)
        //

        CHECK(fixed.isApprox(expected_fixed));
    }

    SECTION("Error Scene 1: A required file (phenotype) is missing")
    {
        // Description: 此测试验证当一个必须的文件（如表型文件）不存在时，
        // 程序会抛出 std::runtime_error 异常以阻止后续执行，确保程序的稳健性。

        // 创建一个有效的 fam 文件，但确保表型文件不存在。
        std::ofstream fam_file(fam_path_);
        fam_file << "FID1 IID1 0 0 1 1\n";
        fam_file.close();

        if (std::filesystem::exists(pheno_path_))
        {
            std::filesystem::remove(pheno_path_);
        }

        REQUIRE_THROWS_AS(
            gelex::DataReader::Create(pheno_path_, fam_path_, "", "", 3, true),
            std::runtime_error);
    }

    SECTION("Error Scene 2: An input file has an invalid header format")
    {
        // Description: 此测试确保文件解析器能正确验证文件头的格式。
        // 如果文件头的前两列不是所要求的 'FID' 和 'IID'，
        // 则应抛出带有明确错误信息的 std::runtime_error 异常。
        create_invalid_header_files();

        REQUIRE_THROWS_WITH(
            gelex::DataReader::Create(pheno_path_, fam_path_, "", "", 3, true),
            Catch::Matchers::ContainsSubstring(
                "First two columns must be 'FID' and 'IID'"));
    }
}
