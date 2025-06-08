#include <string>
#include <vector>

#include <armadillo>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>

#include "gelex/estimator/mcmc_samples.h"
#include "gelex/estimator/mcmc_params.h"
#include "gelex/model/bayes.h"
#include "gelex/model/bayes_effects.h"

using namespace arma;
using namespace gelex;

TEST_CASE("MCMCSamples stores correctly", "[mcmc_samples]")
{
    // Create minimal BayesModel
    dvec phenotype = {1.0, 2.0};
    BayesModel model("y ~ 1 + x", std::move(phenotype));
    
    // Add fixed effect
    std::vector<std::string> names = {"x"};
    std::vector<std::string> levels = {};
    dmat fixed_design = {{0.5}, {1.0}};
    model.add_fixed_effect(std::move(names), std::move(levels), std::move(fixed_design));
    
    // Add random effect
    dmat random_design = eye(2, 2);
    model.add_random_effect("rand", std::move(random_design));
    
    // Add genetic effect
    dmat genetic_design = ones(2, 3);
    model.add_genetic_effect("gen", std::move(genetic_design), BayesAlphabet::A);
    
    // Configure MCMC parameters
    MCMCParams params;
    params.iter = 2000;
    params.n_burnin = 1000;
    params.n_thin = 500;
    params.n_chains = 1;  // Only 1 chain for this test
    
    // Create samples container
    MCMCSamples samples(params, model, params.n_chains);
    
    // Create status with test values
    BayesStatus status(model);
    status.mu.value = 1.0;
    status.fixed.coeff = {0.5};
    status.random[0].coeff = {0.1, 0.2};
    status.random[0].sigma = {0.5};
    status.genetic[0].coeff = {0.1, 0.2, 0.3};
    status.genetic[0].sigma = {0.01};
    status.residual.value = 0.1;
    
    // Store first sample
    samples.store(status, 0, 0);
    
    // Update values and store second sample
    status.mu.value = 2.0;
    status.fixed.coeff = {1.0};
    status.random[0].coeff = {0.3, 0.4};
    status.random[0].sigma = {0.6};
    status.genetic[0].coeff = {0.4, 0.5, 0.6};
    status.genetic[0].sigma = {0.02};
    status.residual.value = 0.2;
    samples.store(status, 1, 0);
    
    // Verify stored values
    // Check dimensions first
    REQUIRE(samples.mu().n_rows == 2);
    REQUIRE(samples.mu().n_cols == 1);
    
    // mu values
    REQUIRE(samples.mu()(0, 0) == 1.0);
    REQUIRE(samples.mu()(1, 0) == 2.0);
    
    // fixed effects
    REQUIRE(samples.fixed().n_rows == 1);  // 1 fixed effect
    REQUIRE(samples.fixed().n_cols == 2);  // 2 samples
    REQUIRE(samples.fixed().n_slices == 1); // 1 chain
    REQUIRE(samples.fixed()(0, 0, 0) == 0.5);
    REQUIRE(samples.fixed()(0, 1, 0) == 1.0);
    
    // random effects
    REQUIRE(samples.random().coeffs[0](0, 0, 0) == 0.1);
    REQUIRE(samples.random().coeffs[0](1, 0, 0) == 0.2);
    REQUIRE(samples.random().sigmas[0](0, 0, 0) == 0.5);
    
    REQUIRE(samples.random().coeffs[0](0, 1, 0) == 0.3);
    REQUIRE(samples.random().coeffs[0](1, 1, 0) == 0.4);
    REQUIRE(samples.random().sigmas[0](0, 1, 0) == 0.6);
    
    // genetic effects
    REQUIRE(samples.genetic().coeffs.size() == 1);  // 1 genetic effect
    REQUIRE(samples.genetic().coeffs[0].n_rows == 3);  // 3 coefficients
    REQUIRE(samples.genetic().coeffs[0].n_cols == 2);  // 2 samples
    REQUIRE(samples.genetic().coeffs[0].n_slices == 1); // 1 chain
    REQUIRE(samples.genetic().coeffs[0](0, 0, 0) == 0.1);
    REQUIRE(samples.genetic().coeffs[0](1, 0, 0) == 0.2);
    REQUIRE(samples.genetic().coeffs[0](2, 0, 0) == 0.3);
    REQUIRE(samples.genetic().sigmas[0](0, 0, 0) == 0.01);
    
    REQUIRE(samples.genetic().coeffs[0](0, 1, 0) == 0.4);
    REQUIRE(samples.genetic().coeffs[0](1, 1, 0) == 0.5);
    REQUIRE(samples.genetic().coeffs[0](2, 1, 0) == 0.6);
    REQUIRE(samples.genetic().sigmas[0](0, 1, 0) == 0.02);
    
    // residual
    REQUIRE(samples.residual()(0, 0) == 0.1);
    REQUIRE(samples.residual()(1, 0) == 0.2);
    
    // Verify h2() calculations
    REQUIRE(samples.h2().n_elem == 2);  // 2 samples
    // h2 = genetic_variance / (genetic_variance + residual_variance)
    double expected_h2_0 = 0.01 / (0.01 + 0.1);
    double expected_h2_1 = 0.02 / (0.02 + 0.2);
    REQUIRE(samples.h2()(0) == Approx(expected_h2_0).margin(1e-10));
    REQUIRE(samples.h2()(1) == Approx(expected_h2_1).margin(1e-10));
}

TEST_CASE("MCMCSamples handles boundary cases", "[mcmc_samples]")
{
    // Create minimal BayesModel
    dvec phenotype = {1.0};
    BayesModel model("y ~ 1", std::move(phenotype));
    
    // Configure MCMC parameters
    MCMCParams params;
    params.iter = 1000;
    params.n_burnin = 1000;
    params.n_thin = 1;
    params.n_chains = 1;
    
    // Create samples container (should have 0 records)
    MCMCSamples samples(params, model, params.n_chains);
    
    // Attempt to store sample (should be ignored)
    BayesStatus status(model);
    status.mu.value = 1.0;
    status.residual.value = 0.1;
    samples.store(status, 0, 0);
    
    // Verify no storage occurred and dimensions are correct
    REQUIRE(samples.mu().n_rows == 0);
    REQUIRE(samples.mu().n_cols == 1);
    REQUIRE(samples.residual().n_rows == 0);
    REQUIRE(samples.residual().n_cols == 1);
    REQUIRE(samples.h2().n_elem == 0);
}
