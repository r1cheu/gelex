#include <string>
#include <vector>

#include <armadillo>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>

#include "gelex/estimator/mcmc_params.h"
#include "gelex/estimator/mcmc_samples.h"
#include "gelex/model/base_effects.h"
#include "gelex/model/bayes.h"
#include "gelex/model/bayes_effects.h"

using namespace arma;   // NOLINT
using namespace gelex;  // NOLINT

TEST_CASE("MCMCSamples stores correctly", "[mcmc_samples]")
{
    // Create minimal BayesModel
    dvec phenotype = {1.0, 2.0};
    BayesModel model("y ~ 1 + x", std::move(phenotype));

    // Add fixed effect
    std::vector<std::string> names = {"x"};
    std::vector<std::string> levels = {"x"};
    dmat fixed_design(2, 1, arma::fill::zeros);
    model.add_fixed_effect(
        std::move(names), std::move(levels), std::move(fixed_design));

    // Add random effect
    dmat random_design = eye(2, 2);
    model.add_random_effect("rand", std::move(random_design));

    // Add genetic effect
    dmat genetic_design = ones(2, 3);
    model.add_genetic_effect(
        "gen", std::move(genetic_design), BayesAlphabet::RR);

    MCMCParams params{};
    params.iter = 2000;
    params.n_burnin = 1000;
    params.n_thin = 500;
    params.n_chains = 2;

    MCMCSamples samples(params, model);

    BayesStatus status(model);
    status.mu.value = 1.0;
    status.fixed.coeff = {0.5};
    status.random[0].coeff = {0.1, 0.2};
    status.random[0].sigma = {0.5};
    status.genetic[0].coeff = {0.1, 0.2, 0.3};
    status.genetic[0].sigma = {0.01};
    status.residual.value = 0.1;

    BayesStatus status_back = status;

    status.mu.value = 2.0;
    status.fixed.coeff = {1.0};
    status.random[0].coeff = {0.3, 0.4};
    status.random[0].sigma = {0.6};
    status.genetic[0].coeff = {0.4, 0.5, 0.6};
    status.genetic[0].sigma = {0.02};
    status.residual.value = 0.2;

    samples.store(status_back, 0, 0);
    samples.store(status, 1, 0);
    samples.store(status, 0, 1);
    samples.store(status_back, 1, 1);

    REQUIRE(samples.mu().n_rows == 2);
    REQUIRE(samples.mu().n_cols == 2);
    REQUIRE(approx_equal(
        samples.mu(), dmat{{1.0, 2.0}, {2.0, 1.0}}, "absdiff", 1e-5));

    dcube fixed(1, 2, 2);
    fixed.slice(0) = rowvec{0.5, 1};
    fixed.slice(1) = rowvec{1, 0.5};
    REQUIRE(samples.fixed().n_rows == 1);    // 1 fixed effect
    REQUIRE(samples.fixed().n_cols == 2);    // 2 samples
    REQUIRE(samples.fixed().n_slices == 2);  // 2 chain
    REQUIRE(approx_equal(samples.fixed(), fixed, "absdiff", 1e-5));

    dcube random(2, 2, 2);
    random.slice(0) = dmat{{0.1, 0.3}, {0.2, 0.4}};
    random.slice(1) = dmat{{0.3, 0.1}, {0.4, 0.2}};

    REQUIRE(samples.random().coeffs.size() == 1);  // 2 random effects
    REQUIRE(samples.random().coeffs[0].n_rows == 2);
    REQUIRE(samples.random().coeffs[0].n_cols == 2);
    REQUIRE(samples.random().coeffs[0].n_slices == 2);
    REQUIRE(approx_equal(samples.random().coeffs[0], random, "absdiff", 1e-5));

    dcube random_sigma(1, 2, 2);
    random_sigma.slice(0) = rowvec{0.5, 0.6};
    random_sigma.slice(1) = rowvec{0.6, 0.5};
    REQUIRE(samples.random().sigmas.size() == 1);  // 1 random effect
    REQUIRE(samples.random().sigmas[0].n_rows == 1);
    REQUIRE(samples.random().sigmas[0].n_cols == 2);
    REQUIRE(samples.random().sigmas[0].n_slices == 2);
    REQUIRE(approx_equal(
        samples.random().sigmas[0], random_sigma, "absdiff", 1e-5));

    dcube genetic(3, 2, 2);
    REQUIRE(samples.genetic().coeffs.size() == 1);  // 1 genetic effect
    REQUIRE(samples.genetic().coeffs[0].n_rows == 3);
    REQUIRE(samples.genetic().coeffs[0].n_cols == 2);
    REQUIRE(samples.genetic().coeffs[0].n_slices == 2);

    genetic.slice(0) = dmat{{0.1, 0.4}, {0.2, 0.5}, {0.3, 0.6}};
    genetic.slice(1) = dmat{{0.4, 0.1}, {0.5, 0.2}, {0.6, 0.3}};
    REQUIRE(
        approx_equal(samples.genetic().coeffs[0], genetic, "absdiff", 1e-5));

    dcube genetic_sigma(1, 2, 2);
    REQUIRE(samples.genetic().sigmas.size() == 1);  // 1 genetic effect
    REQUIRE(samples.genetic().sigmas[0].n_rows == 1);
    REQUIRE(samples.genetic().sigmas[0].n_cols == 2);
    REQUIRE(samples.genetic().sigmas[0].n_slices == 2);

    genetic_sigma.slice(0) = rowvec{0.01, 0.02};
    genetic_sigma.slice(1) = rowvec{0.02, 0.01};
    REQUIRE(approx_equal(
        samples.genetic().sigmas[0], genetic_sigma, "absdiff", 1e-5));

    dmat residual(1, 2);

    REQUIRE(samples.residual().n_rows == 2);
    REQUIRE(samples.residual().n_cols == 2);
    REQUIRE(approx_equal(
        samples.residual(), dmat{{0.1, 0.2}, {0.2, 0.1}}, "absdiff", 1e-5));
}

TEST_CASE("MCMCSamples handles boundary cases", "[mcmc_samples]")
{
    // Create minimal BayesModel
    dvec phenotype = {1.0};
    BayesModel model("y ~ 1", std::move(phenotype));

    // Configure MCMC parameters
    MCMCParams params{};
    params.iter = 1000;
    params.n_burnin = 1000;
    params.n_thin = 1;
    params.n_chains = 1;

    // Create samples container (should have 0 records)
    MCMCSamples samples(params, model);

    // Attempt to store sample (should be ignored)
    BayesStatus status(model);
    status.mu.value = 1.0;
    status.residual.value = 0.1;
    samples.store(status, 0, 0);

    // Verify no storage occurred and dimensions are correct
    REQUIRE(samples.mu().n_rows == 0);
    REQUIRE(samples.mu().n_cols == 1);
    REQUIRE(samples.fixed().empty());
    REQUIRE(samples.random().coeffs.empty());
    REQUIRE(samples.random().sigmas.empty());
    REQUIRE(samples.genetic().coeffs.empty());
    REQUIRE(samples.genetic().sigmas.empty());
    REQUIRE(samples.residual().n_rows == 0);
    REQUIRE(samples.residual().n_cols == 1);
    REQUIRE(samples.h2().n_elem == 0);
}
