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

#include <cmath>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/algo/infer/params.h"
#include "gelex/algo/infer/vi/solver.h"
#include "gelex/algo/infer/vi/trait_model.h"
#include "gelex/data/genotype/genotype_matrix.h"
#include "gelex/model/bayes/genotype_storage.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_config.h"
#include "gelex/types/bayes_method.h"
#include "gelex/types/fixed_effects.h"
#include "gelex/types/vi_result.h"

using namespace gelex;  // NOLINT

namespace
{

auto make_geno_matrix(Eigen::Index n_samples, Eigen::Index n_snps)
    -> GenotypeMatrix
{
    Eigen::MatrixXd data = Eigen::MatrixXd::Random(n_samples, n_snps);
    auto mean = data.colwise().mean().transpose().eval();
    auto stddev = Eigen::VectorXd::Ones(n_snps);
    return GenotypeMatrix(
        std::move(data), {}, std::move(mean), std::move(stddev));
}

auto make_bayes_rr_model(Eigen::Index n_samples, Eigen::Index n_snps)
    -> std::pair<BayesModel, bayes::Priors>
{
    auto phenotype = Eigen::VectorXd::Random(n_samples);
    auto fixed = FixedEffect::build(n_samples);
    auto geno = make_geno_matrix(n_samples, n_snps);

    std::vector<bayes::GeneticEffect> genetics;
    genetics.emplace_back(
        GeneticMode::A, bayes::GenotypeStorage{std::move(geno)});

    double pheno_var
        = phenotype.array().square().mean() - std::pow(phenotype.mean(), 2.0);
    PriorSetConfig pc(BayesMethodConfig{BayesBase::RR}, pheno_var);
    bayes::Priors priors(pc, genetics, 0);

    BayesModel model(phenotype, std::move(fixed), std::move(genetics));
    return {std::move(model), std::move(priors)};
}

}  // namespace

TEST_CASE("VIState construction from BayesRR model", "[cavi]")
{
    constexpr Eigen::Index kN = 20;
    constexpr Eigen::Index kP = 5;
    auto [model, priors] = make_bayes_rr_model(kN, kP);

    vi::State state(model, priors);

    // Fixed: intercept only
    CHECK(state.fixed().coeffs.size() == 1);

    // One genetic effect (additive)
    REQUIRE(state.genetics().size() == 1);
    const auto& gs = state.genetics()[0];
    CHECK(gs.type == GeneticMode::A);
    CHECK(gs.coeffs.size() == kP);
    CHECK(gs.sigma2.size() == kP);
    CHECK(gs.u.size() == kN);
    CHECK(gs.marker_variance.size() == kP);

    // Residual
    CHECK(state.residual().y_adj.size() == kN);
    CHECK(state.residual().variance > 0.0);
}

TEST_CASE("CAVI RR single iteration produces correct posteriors", "[cavi]")
{
    constexpr Eigen::Index kN = 10;
    constexpr Eigen::Index kP = 3;

    // Deterministic genotype matrix
    Eigen::MatrixXd X(kN, kP);
    X << 1, 0, 2, 0, 1, 1, 2, 2, 0, 1, 0, 1, 0, 1, 2, 2, 0, 1, 1, 1, 0, 0, 2, 1,
        1, 1, 0, 2, 0, 2;

    Eigen::VectorXd mean = X.colwise().mean().transpose();
    auto stddev = Eigen::VectorXd::Ones(kP);
    GenotypeMatrix geno(
        Eigen::MatrixXd(X), {}, std::move(mean), std::move(stddev));

    Eigen::VectorXd phenotype = Eigen::VectorXd::LinSpaced(kN, -1.0, 1.0);

    auto fixed = FixedEffect::build(kN);
    std::vector<bayes::GeneticEffect> genetics;
    genetics.emplace_back(
        GeneticMode::A, bayes::GenotypeStorage{std::move(geno)});

    double pheno_var
        = phenotype.array().square().mean() - std::pow(phenotype.mean(), 2.0);
    PriorSetConfig pc(BayesMethodConfig{BayesBase::RR}, pheno_var);
    bayes::Priors priors(pc, genetics, 0);

    BayesModel model(phenotype, std::move(fixed), std::move(genetics));
    vi::State state(model, priors);

    const auto& effect = model.genetics()[0];
    const auto* prior = priors.genetic(GeneticMode::A);
    REQUIRE(prior != nullptr);

    const auto& csn = effect.XtX_diag;
    const double res_var = state.residual().variance;
    const double marker_var = state.genetics()[0].marker_variance(0);

    // Run one iteration of the RR updater
    vi::RR trait_model{};
    trait_model(model, priors, state);

    const auto& gs = state.genetics()[0];

    // After one iteration, for each marker sequentially updated:
    //   post_var_j = σ²_e / (X_j'X_j + σ²_e / σ²_β)
    //   post_mean_j = (X_j'·y_adj + X_j'X_j·old_j) / (X_j'X_j + σ²_e / σ²_β)
    // Since markers update sequentially (y_adj changes), we can verify
    // that post_var = sigma2 for each marker
    for (Eigen::Index j = 0; j < kP; ++j)
    {
        const double v = csn(j) + res_var / marker_var;
        const double expected_var = res_var / v;
        CHECK(gs.sigma2(j) > 0.0);
        CHECK_THAT(
            gs.sigma2(j), Catch::Matchers::WithinRel(expected_var, 1e-10));
    }

    // Breeding values should be updated
    Eigen::VectorXd expected_u = bayes::get_matrix_ref(effect.X) * gs.coeffs;
    CHECK(gs.u.isApprox(expected_u, 1e-10));
}

TEST_CASE("CAVI RR converges on synthetic data", "[cavi]")
{
    constexpr Eigen::Index kN = 30;
    constexpr Eigen::Index kP = 3;

    // Deterministic design: LinSpaced columns, orthogonal-ish
    Eigen::MatrixXd geno_data(kN, kP);
    geno_data.col(0) = Eigen::VectorXd::LinSpaced(kN, -1.0, 1.0);
    geno_data.col(1) = Eigen::VectorXd::LinSpaced(kN, 0.5, -0.5);
    geno_data.col(2) = Eigen::VectorXd::LinSpaced(kN, -0.3, 0.3);

    Eigen::VectorXd mean = geno_data.colwise().mean().transpose();
    auto stddev = Eigen::VectorXd::Ones(kP);
    GenotypeMatrix geno(
        Eigen::MatrixXd(geno_data), {}, std::move(mean), std::move(stddev));

    Eigen::VectorXd beta_true(kP);
    beta_true << 1.0, -0.5, 0.3;

    Eigen::VectorXd phenotype = geno_data * beta_true;

    auto fixed = FixedEffect::build(kN);
    std::vector<bayes::GeneticEffect> genetics;
    genetics.emplace_back(
        GeneticMode::A, bayes::GenotypeStorage{std::move(geno)});

    double pheno_var
        = phenotype.array().square().mean() - std::pow(phenotype.mean(), 2.0);
    PriorSetConfig pc(BayesMethodConfig{BayesBase::RR}, pheno_var);
    bayes::Priors priors(pc, genetics, 0);

    BayesModel model(phenotype, std::move(fixed), std::move(genetics));

    // Collect progress events
    std::vector<double> deltas;
    auto observer = FitObserver(
        [&](const FitEvent& event)
        {
            if (const auto* p = std::get_if<FitVIProgressEvent>(&event))
            {
                if (!p->done)
                {
                    deltas.push_back(p->delta);
                }
            }
        });

    vi::Params params;
    params.max_iters = 200;
    params.tol = 1e-8;

    vi::Solver cavi(params, vi::RR{});
    auto result = cavi.run(model, priors, observer);

    // Should have converged before max_iters
    CHECK(deltas.size() < 200);

    // Final delta should be below tolerance
    CHECK(deltas.back() < params.tol);

    // Estimated coefficients should correlate with true effects
    const auto* gs = result.genetic(GeneticMode::A);
    REQUIRE(gs != nullptr);

    // Ridge regression shrinks coefficients toward zero, so exact match is
    // not expected. Verify the direction is correct (positive correlation).
    const double corr = (gs->coeffs.array() * beta_true.array()).sum()
                        / (gs->coeffs.norm() * beta_true.norm());
    CHECK(corr > 0.8);
}
