#include "chenx/predictor.h"

#include <memory>

#include <armadillo>

#include "chenx/data/cross_grm.h"
#include "chenx/model/linear_mixed_model.h"

namespace chenx
{
using arma::rowvec;
void Predictor::set_cross_grm(
    std::string_view method,
    rowvec&& center,
    double scale_factor,
    uint64_t chunk_size)
{
    if (method == "add")
    {
        cross_grms_.emplace_back(std::make_unique<AddCrossGrm>(
            train_bed_,
            std::move(center),
            scale_factor,
            chunk_size,
            params_.dropped_individuals()));
    }
    else if (method == "dom")
    {
        cross_grms_.emplace_back(std::make_unique<DomCrossGrm>(
            train_bed_,
            std::move(center),
            scale_factor,
            chunk_size,
            params_.dropped_individuals()));
    }
    else
    {
        throw std::runtime_error("Unknown method: " + std::string(method));
    }
};

void Predictor::set_grm(dmat&& grm)
{
    grms_.emplace_back(std::move(grm));
}

std::pair<dmat, dvec>
Predictor::solver_chol(dmat& V, const dmat& X, const dvec& y)
{
    char uplo = 'L';
    int n = static_cast<int>(V.n_cols);
    int info{};
    arma::lapack::potrf(&uplo, &n, V.memptr(), &n, &info);
    if (info != 0)
    {
        throw std::runtime_error("V Matrix is not symmetric positive definite");
    }

    dmat rhs = arma::join_horiz(X, y);
    int m = static_cast<int>(rhs.n_cols);
    arma::lapack::potrs(&uplo, &n, &m, V.memptr(), &n, rhs.memptr(), &n, &info);

    dmat phi_1 = rhs.cols(0, X.n_cols - 1);
    dvec phi_2 = rhs.col(X.n_cols);

    return std::make_pair(phi_1, phi_2);
}

dvec Predictor::ComputePy()
{
    dmat V{grms_[0] * params_.sigma()[0]};
    for (size_t i = 1; i < grms_.size(); ++i)
    {
        V += grms_[i] * params_.sigma()[i];
    }
    auto [phi_1, phi_2] = solver_chol(V, params_.X(), params_.y());
    dmat X_t{params_.X().t()};
    return phi_2 - phi_1 * arma::inv(X_t * phi_1) * X_t * phi_2;
}

dmat Predictor::predict(std::string_view test_bed)
{
    return cross_grms_[0]->Compute(test_bed);
};

/*dmat Predictor::predict(std::string_view test_bed)*/
/*{*/
/*    if (Py_.empty())*/
/*    {*/
/*        Py_ = ComputePy();*/
/*    }*/
/*    dmat U;*/
/*    for (int i = 0; i < cross_grms_.size(); ++i)*/
/*    {*/
/*        dmat new_k = cross_grms_[i]->Compute(test_bed);*/
/*        dvec u = new_k * Py_ * params_.sigma()[i];*/
/*        if (U.empty())*/
/*        {*/
/*            U = u;*/
/*        }*/
/*        else*/
/*        {*/
/*            U = arma::join_horiz(U, u);*/
/*        }*/
/*    }*/
/*    return U;*/
/*};*/
}  // namespace chenx
