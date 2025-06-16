#pragma once
#include <memory>
#include <string_view>
#include <vector>

#include <armadillo>

#include "gelex/data/cross_grm.h"
#include "gelex/model/freq/model.h"

namespace gelex
{
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
        arma::rowvec center,
        double scale_factor,
        size_t chunk_size);
    arma::dmat compute_random_effects(std::string_view test_bed);
    arma::dmat compute_fixed_effects(
        const arma::dvec& covariates) const noexcept;

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
