#pragma once
#include <string_view>
#include <vector>

#include <armadillo>

#include "gelex/data/grm.h"
#include "gelex/model/freq/model.h"

namespace gelex
{

struct BLUP
{
    arma::dvec level_solutions;  // Renamed from u to level_solutions
    std::string name;
};

class GBLUPPredictor
{
   public:
    GBLUPPredictor(
        std::string_view train_bed,
        const std::vector<std::string>& train_id_order,
        const GBLUP& model);
    GBLUPPredictor(
        std::string_view train_bed,
        const std::vector<std::string>& train_id_order,
        const arma::dvec& beta,
        const std::vector<BLUP>& blups);

    GBLUPPredictor(const GBLUPPredictor&) = delete;
    GBLUPPredictor(GBLUPPredictor&&) noexcept = default;
    const GBLUPPredictor& operator=(const GBLUPPredictor&) = delete;
    GBLUPPredictor& operator=(GBLUPPredictor&&) noexcept = default;

    ~GBLUPPredictor() = default;

    void add_cross_grm(
        std::string_view method,
        arma::dvec p_major,
        double scale_factor,
        size_t chunk_size);

    arma::dmat compute_random_effects(std::string_view test_bed);
    arma::dmat compute_fixed_effects(const arma::dvec& covariates) const;

   private:
    std::string train_bed_;
    std::vector<std::string> train_id_order_;
    std::vector<CrossGRM> cross_grms_;

    arma::dvec beta_;
    std::vector<BLUP> blups_;
};
}  // namespace gelex
