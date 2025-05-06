#include <string>

#include <armadillo>
#include <catch2/catch_test_macros.hpp>
#include <utility>

#include "gelex/estimator/estimator.h"
#include "gelex/model/gblup.h"

TEST_CASE("Linear Mixed Model Fitted Check")
{
    using arma::dmat;
    using arma::dvec;

    arma::dvec Phenotype;
    Phenotype.load(
        std::string(GELEX_TESTS_DIR) + "/data/wheat100_phenotype.bin");
    arma::dmat A;
    A.load(std::string(GELEX_TESTS_DIR) + "/data/wheat100_grm.bin");
    arma::dmat X = arma::ones<arma::dmat>(Phenotype.n_elem, 1);
    gelex::GBLUP model{"yield~1", std::move(Phenotype)};
    model.add_fixed_effect(
        std::vector<std::string>{"intercept"},
        std::vector<std::string>{"intercept"},
        std::move(X));
    model.add_genetic_effect("A", arma::speye(arma::size(A)), std::move(A));

    dvec sigma_hat{0.161, 0.513};

    SECTION("Average Information")
    {
        gelex::Estimator estimator{"AI", 20, 1e-6};
        estimator.fit(model, false, true);
        REQUIRE(arma::approx_equal(model.sigma(), sigma_hat, "absdiff", 1e-3));
    }
}
