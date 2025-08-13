#pragma once
#include <string_view>
#include <vector>

#include <armadillo>

#include "gelex/data/grm.h"
#include "gelex/model/freq/model.h"

namespace gelex
{

class GBLUPPredictor
{
   public:
    GBLUPPredictor(
        std::string_view train_bed,
        const std::vector<std::string>& train_id_order,
        const GBLUP& model);

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

    arma::dmat compute_fixed_effects(const arma::dvec& covariates) const;
    // arma::dmat compute_genetic_effects(std::string_view test_bed);
    // arma::dvec compute_random_effects(
    //     std::string_view name,
    //     const arma::dmat& design_matrix) const;

   private:
    struct Effect
    {
        std::string name;
        arma::sp_dmat design_matrix;
        double sigma;
    };
    std::string train_bed_;
    std::vector<std::string> train_id_order_;
    std::vector<CrossGRM> cross_grms_;

    arma::dvec beta_;
    std::vector<Effect> blups_;
};
}  // namespace gelex
