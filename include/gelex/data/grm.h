#pragma once

#include <armadillo>
#include <cstddef>

#include "gelex/data/bed_reader.h"

namespace gelex
{

class GRM
{
   public:
    explicit GRM(
        std::string_view bed_file,
        size_t chunk_size = 10000,
        const std::vector<std::string>& target_order = {});
    GRM(const GRM&) = delete;
    GRM(GRM&&) noexcept = default;
    GRM& operator=(const GRM&) = delete;
    GRM& operator=(GRM&&) noexcept = default;

    ~GRM() = default;

    arma::dmat compute(bool add = true);
    arma::dvec p_major() const noexcept { return p_major_; }
    double scale_factor() const { return scale_factor_; }

    const std::vector<std::string>& individuals() const noexcept
    {
        return bed_.individuals();
    }

   private:
    BedReader bed_;
    double scale_factor_{1.0};
    arma::dvec p_major_;
};

class CrossGRM
{
   public:
    explicit CrossGRM(
        std::string_view train_bed,
        arma::dvec p_major,
        double scale_factor,
        size_t chunk_size = 10000,
        const std::vector<std::string>& target_order = {});

    CrossGRM(const CrossGRM&) = delete;
    CrossGRM(CrossGRM&&) noexcept = default;
    CrossGRM& operator=(const CrossGRM&) = delete;
    CrossGRM& operator=(CrossGRM&&) noexcept = default;

    ~CrossGRM() = default;

    const std::vector<std::string>& individuals() const noexcept
    {
        return individuals_;
    }

    arma::dmat compute(std::string_view test_bed, bool add = true);

   private:
    BedReader bed_;
    void check_snp_consistency(const BedReader& test_bed) const;
    arma::dvec p_major_;
    double scale_factor_;

    std::vector<std::string> individuals_;
    size_t chunk_size_;
};

arma::dvec compute_p_major(const arma::dmat& genotype);
void code_method_varden(
    arma::dvec p_major,
    arma::dmat& genotype,
    bool add = true);

}  // namespace gelex
