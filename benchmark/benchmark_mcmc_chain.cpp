/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <random>
#include <vector>

#include <nanobench.h>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/algo/infer/mcmc/recipes.h"
#include "gelex/data/genotype/matrix.h"
#include "gelex/data/genotype/processor.h"
#include "gelex/data/genotype/storage.h"
#include "gelex/model/bayes/builder.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/method.h"
#include "gelex/model/bayes/model.h"
#include "gelex/types/fixed_effects.h"
#include "gelex/types/genetic_effect_type.h"

namespace
{

using gelex::BayesBase;
using gelex::BayesModel;
using gelex::FixedEffect;
using gelex::GeneticMode;
using gelex::PriorOverrides;
using gelex::bayes::BayesConfig;
using gelex::bayes::DominancePolicy;
using gelex::bayes::GeneticEffect;
using gelex::bayes::GenotypeStorage;
using gelex::genotype::GenotypeMatrix;

constexpr Eigen::Index kIndividuals = 100;
constexpr Eigen::Index kMarkers = 200;

auto make_genotype(uint64_t seed) -> GenotypeMatrix
{
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uniform(0.05, 0.5);

    Eigen::MatrixXd X(kIndividuals, kMarkers);
    for (Eigen::Index j = 0; j < kMarkers; ++j)
    {
        std::binomial_distribution<int> dist(2, uniform(rng));
        X.col(j) = Eigen::VectorXd::NullaryExpr(
            kIndividuals,
            [&](Eigen::Index) { return static_cast<double>(dist(rng)); });
    }

    Eigen::VectorXd mean(kMarkers);
    Eigen::VectorXd stddev(kMarkers);
    std::vector<int64_t> mono;
    auto process = gelex::genotype::Standardize<GeneticMode::A>::process;
    for (Eigen::Index j = 0; j < kMarkers; ++j)
    {
        auto col = X.col(j);
        const auto stats = process(col);
        mean(j) = stats.mean;
        stddev(j) = stats.stddev;
        if (stats.is_monomorphic)
        {
            mono.push_back(static_cast<int64_t>(j));
        }
    }

    return GenotypeMatrix(
        std::move(X), std::move(mono), std::move(mean), std::move(stddev));
}

auto make_phenotype(const Eigen::MatrixXd& X, uint64_t seed) -> Eigen::VectorXd
{
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    Eigen::VectorXd betas = Eigen::VectorXd::Zero(X.cols());
    std::uniform_int_distribution<Eigen::Index> idx(0, X.cols() - 1);
    for (int k = 0; k < 30; ++k)
    {
        betas(idx(rng)) = normal(rng) * 0.25;
    }
    return (X * betas).unaryExpr([&](double v) { return v + normal(rng); });
}

auto make_model() -> BayesModel
{
    auto geno = make_genotype(0xBEEF1234ULL);
    auto y = make_phenotype(geno.matrix(), 0xCAFE5678ULL);
    auto fixed = FixedEffect::build(kIndividuals);
    std::vector<GeneticEffect> genetics;
    genetics.emplace_back(GeneticMode::A, GenotypeStorage{std::move(geno)});
    return BayesModel(std::move(y), std::move(fixed), std::move(genetics));
}

template <typename Factory>
void bench_chain(
    ankerl::nanobench::Bench& b,
    const BayesModel& model,
    const char* name,
    BayesConfig cfg,
    Factory factory,
    Eigen::Index n_iters)
{
    auto method = gelex::build_bayes_method(
        PriorOverrides{
            .method = cfg,
            .phenotype_variance = model.phenotype_variance(),
        },
        model);
    gelex::mcmc::State state{model, method};
    std::mt19937_64 rng(42);
    gelex::mcmc::Context ctx{model, method, state, rng};
    auto chain = factory(ctx);
    b.run(
        name,
        [&]
        {
            for (Eigen::Index i = 0; i < n_iters; ++i)
            {
                chain.step();
            }
            ankerl::nanobench::doNotOptimizeAway(state.residual().variance);
        });
}

}  // namespace

TEST_CASE("MCMC trait chain per-iteration", "[!benchmark][mcmc][chain]")
{
    const auto model = make_model();
    constexpr Eigen::Index kIters = 50;

    ankerl::nanobench::Bench b;
    b.title(
         "MCMC trait per-iter (n=" + std::to_string(kIndividuals)
         + ", m=" + std::to_string(kMarkers) + ")")
        .unit("iter")
        .batch(static_cast<double>(kIters))
        .warmup(1)
        .minEpochIterations(20);

    using gelex::GeneticMode;
    constexpr auto kSym = DominancePolicy::symmetric;
    const BayesConfig cfg_a{BayesBase::A, GeneticMode::A, kSym, false};
    const BayesConfig cfg_b{BayesBase::B, GeneticMode::A, kSym, false};
    const BayesConfig cfg_c{BayesBase::C, GeneticMode::A, kSym, false};
    const BayesConfig cfg_rr{BayesBase::RR, GeneticMode::A, kSym, false};
    const BayesConfig cfg_bpi{BayesBase::B, GeneticMode::A, kSym, true};
    const BayesConfig cfg_cpi{BayesBase::C, GeneticMode::A, kSym, true};

    bench_chain(
        b,
        model,
        "A (chain)",
        cfg_a,
        gelex::mcmc::make_bayes_a_chain<GeneticMode::A>,
        kIters);
    bench_chain(
        b,
        model,
        "B (chain)",
        cfg_b,
        gelex::mcmc::make_bayes_b_chain<GeneticMode::A>,
        kIters);
    bench_chain(
        b,
        model,
        "C (chain)",
        cfg_c,
        gelex::mcmc::make_bayes_c_chain<GeneticMode::A>,
        kIters);
    bench_chain(
        b,
        model,
        "RR (chain)",
        cfg_rr,
        gelex::mcmc::make_bayes_rr_chain<GeneticMode::A>,
        kIters);
    bench_chain(
        b,
        model,
        "Bpi (chain)",
        cfg_bpi,
        gelex::mcmc::make_bayes_bpi_chain<GeneticMode::A>,
        kIters);
    bench_chain(
        b,
        model,
        "Cpi (chain)",
        cfg_cpi,
        gelex::mcmc::make_bayes_cpi_chain<GeneticMode::A>,
        kIters);
}
