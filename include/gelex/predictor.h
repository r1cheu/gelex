#pragma once
#include <memory>
#include <string_view>
#include <vector>

#include <armadillo>

#include "gelex/data/cross_grm.h"
#include "gelex/model/gblup.h"

namespace gelex
{
using arma::rowvec;
class Predictor
{
   public:
    Predictor(std::string_view train_bed, GBLUPParams params)
        : train_bed_{train_bed}, params_{std::move(params)} {};
    Predictor(const Predictor&) = delete;
    Predictor(Predictor&&) noexcept = default;
    const Predictor& operator=(const Predictor&) = delete;
    Predictor& operator=(Predictor&&) noexcept = default;

    ~Predictor() = default;

    void set_cross_grm(
        std::string_view method,
        rowvec center,
        double scale_factor,
        uint64_t chunk_size);
    dmat ComputeRandomEffects(std::string_view test_bed);
    dmat ComputeFixedEffects(const dvec& covariates) const noexcept;

    const std::vector<std::string>& test_individuals() const noexcept
    {
        return cross_grms_.front()->test_individuals();
    }

   private:
    std::string train_bed_;
    std::vector<std::unique_ptr<CrossGrm>> cross_grms_;
    GBLUPParams params_;
};
}  // namespace gelex
