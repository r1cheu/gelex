#include <armadillo>

#include "chenx/estimator.h"
#include "chenx/model/linear_mixed_model.h"

int main(int argc, char* argv[])
{
    arma::dvec y;
    y.load("/home/rlchen/project/phenx/tests/wheat100_phenotype.bin");
    arma::dcube A(100, 100, 1);
    A.slice(0).load("/home/rlchen/project/phenx/tests/wheat100_grm.bin");
    arma::dmat X(100, 1, arma::fill::ones);

    // Create linear mixed model
    chenx::LinearMixedModel model{
        std::move(y),
        std::move(X),
        std::move(A),
        std::vector<std::string>{"a"}};

    chenx::Estimator estimator{"NR", 20, 1e-6};

    estimator.Fit(model);
    estimator.set_optimizer("FS", 20, 1e-6);
    estimator.Fit(model);
    estimator.set_optimizer("AI", 20, 1e-6);
    estimator.Fit(model);
    return 0;
}
