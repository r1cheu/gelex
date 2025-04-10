#include <string>

#include <armadillo>
#include <catch2/catch_test_macros.hpp>

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
    arma::dcube rands{Phenotype.n_elem, Phenotype.n_elem, 1};
    rands.slice(0) = A;

    gelex::GBLUP model{
        std::move(Phenotype),
        std::move(X),
        std::move(rands),
        std::vector<std::string>{"random"}};

    dvec sigma_hat{0.161, 0.513};

    SECTION("Newton Raphson")
    {
        gelex::Estimator estimator{"NR", 20, 1e-6};
        estimator.Fit(model, true, false);
        REQUIRE(arma::approx_equal(model.sigma(), sigma_hat, "absdiff", 1e-3));
    }

    SECTION("Fisher Scoring")
    {
        gelex::Estimator estimator{"FS", 20, 1e-6};
        estimator.Fit(model, true, false);
        REQUIRE(arma::approx_equal(model.sigma(), sigma_hat, "absdiff", 1e-3));
    }

    SECTION("Average Information")
    {
        gelex::Estimator estimator{"AI", 20, 1e-6};
        estimator.Fit(model, true, false);
        REQUIRE(arma::approx_equal(model.sigma(), sigma_hat, "absdiff", 1e-3));
    }
}
