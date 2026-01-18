#include <benchmark/benchmark.h>
#include <cstdio>
#include "gelex/gwas/gwas_writer.h"
#include "gelex/types/snp_info.h"

using namespace gelex::gwas;

// 创建一个写往 /dev/null 的 fixture，避免磁盘 I/O 波动影响格式化性能测试
class GwasWriterFixture : public benchmark::Fixture
{
   public:
    std::string temp_filename;
    gelex::SnpMeta snp;
    GwasWriter::AssocResult result;

    void SetUp(const benchmark::State&)
    {
// Linux/Mac 使用 /dev/null, Windows 使用 NUL
#ifdef _WIN32
        temp_filename = "NUL";
#else
        temp_filename = "/dev/null";
        // 注意：GwasWriter 会自动添加 .gwas.tsv 后缀，
        // 所以实际路径会变成 /dev/null.gwas.tsv，这在 Linux
        // 下可能需要权限或会报错。 为了测试纯代码性能，建议暂时 hack 一下
        // GwasWriter 或者输出到 /tmp/bench_test
        temp_filename = "/tmp/bench_test_ignore";
#endif

        // 初始化测试数据
        snp = {"1", "rs123456", 100000, 'A', 'G'};
        result.beta = 0.0123;
        result.se = 0.0045;
        result.p_value = 1.23e-8;
    }

    void TearDown(const benchmark::State&)
    {
        // 清理文件
        std::remove((temp_filename + ".gwas.tsv").c_str());
    }
};

BENCHMARK_DEFINE_F(GwasWriterFixture, BM_WriteResult)(benchmark::State& state)
{
    // 使用 AdditiveDominance 模型，因为它写入的字段最多，压力最大
    GwasWriter writer(temp_filename);
    writer.write_header();

    for (auto _ : state)
    {
        writer.write_result(snp, result);
    }
}

BENCHMARK_REGISTER_F(GwasWriterFixture, BM_WriteResult);

BENCHMARK_MAIN();
