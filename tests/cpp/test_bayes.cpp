#include <random>

#include <catch2/catch_test_macros.hpp>

#include "gelex/estimator/bayes/mcmc.h"
#include "gelex/model/bayes/model.h"

struct SimulationData
{
    arma::mat genotypes;     // 基因型矩阵 X (n_samples x n_markers)
    arma::vec phenotypes;    // 表型向量 y (n_samples x 1)
    arma::vec true_effects;  // 真实的遗传效应 beta (n_markers x 1)
};

SimulationData generate_simulation_data(
    const arma::uword n_samples,
    const arma::uword n_markers,
    const arma::uword n_causal,
    const double heritability)
{
    if (n_causal > n_markers)
    {
        throw std::invalid_argument(
            "Causal markers cannot exceed total markers.");
    }

    if (heritability < 0.0 || heritability > 1.0)
    {
        throw std::invalid_argument("Heritability must be between 0 and 1.");
    }

    // 1. 初始化随机数生成器
    std::mt19937_64 gen(42);

    arma::mat genotypes(n_samples, n_markers);
    std::uniform_real_distribution<> maf_dist(0.05, 0.5);

    for (arma::uword j = 0; j < n_markers; ++j)
    {
        double maf = maf_dist(gen);
        std::binomial_distribution<> geno_dist(2, maf);
        for (arma::uword i = 0; i < n_samples; ++i)
        {
            genotypes(i, j) = static_cast<double>(geno_dist(gen));
        }
    }

    // 3. 生成真实的遗传效应向量 (beta)
    // 只有 n_causal 个标记具有非零效应
    arma::vec true_effects = arma::zeros<arma::vec>(n_markers);

    // 随机选择致因标记
    std::vector<arma::uword> marker_indices(n_markers);
    std::iota(marker_indices.begin(), marker_indices.end(), 0);
    std::shuffle(marker_indices.begin(), marker_indices.end(), gen);

    // 从标准正态分布 N(0, 1) 中为致因标记抽取效应值
    std::normal_distribution<> effect_dist(0.0, 1.0);
    for (arma::uword i = 0; i < n_causal; ++i)
    {
        true_effects(marker_indices[i]) = effect_dist(gen);
    }

    // 4. 计算遗传值 (g = X * beta) 并标准化
    arma::vec genetic_value = genotypes * true_effects;

    // 5. 生成表型 (y = g + epsilon)
    // 根据遗传力 h^2 计算遗传方差和环境方差
    // h^2 = Var(g) / Var(y) = Var(g) / (Var(g) + Var(e))
    double var_g = arma::var(genetic_value);
    arma::vec phenotypes;

    if (heritability == 1.0)
    {
        phenotypes = genetic_value;
    }
    else if (heritability == 0.0)
    {
        // 如果遗传力为0，表型完全由噪声决定
        // 为了有一个合理的尺度，我们让噪声的方差等于1
        phenotypes = arma::randn<arma::vec>(n_samples);
    }
    else
    {
        double var_e = var_g * (1.0 - heritability) / heritability;
        double sd_e = std::sqrt(var_e);

        // 生成环境噪声 epsilon ~ N(0, var_e)
        arma::vec epsilon = arma::randn<arma::vec>(n_samples) * sd_e;

        phenotypes = genetic_value + epsilon;
    }

    // 返回所有生成的数据
    return {genotypes, phenotypes, true_effects};
}
using namespace arma;

TEST_CASE("Bayes MCMC Estimation", "[bayes][mcmc]")
{
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> maf_dist(0.05, 0.5);
    dmat X(100, 1000);
    for (arma::uword j = 0; j < 1000; ++j)
    {
        double maf = maf_dist(rng);
        std::binomial_distribution<int> geno_dist(2, maf);
        for (arma::uword i = 0; i < 100; ++i)
        {
            X.at(i, j) = static_cast<double>(geno_dist(rng));
        }
    }

    SECTION("BayesRR") {}
}
