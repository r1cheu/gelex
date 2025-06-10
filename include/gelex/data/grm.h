#pragma once

#include <armadillo>

#include <string_view>
#include <vector>
#include "gelex/data/bed_reader.h"

namespace gelex
{

void dom_encode(arma::dmat& genotype);

class IGrm
{
   public:
    explicit IGrm(
        std::string_view bed_file,
        size_t chunk_size = DEFAULT_CHUNK_SIZE,
        const std::vector<std::string>& exclude_individuals = {});
    IGrm(const IGrm&) = delete;
    IGrm(IGrm&&) noexcept = default;
    IGrm& operator=(const IGrm&) = delete;
    IGrm& operator=(IGrm&&) noexcept = default;
    virtual ~IGrm() = default;

    double scale_factor() const noexcept { return scale_factor_; }
    void set_scale_factor(double scale_factor) { scale_factor_ = scale_factor; }

    const BedReader& bed() const noexcept { return bed_; }
    BedReader& bed() noexcept { return bed_; }

    const arma::rowvec& center() const noexcept { return center_; }
    void set_center(arma::rowvec&& center) { center_ = std::move(center); }
    void set_center(size_t start, const arma::rowvec& center)
    {
        center_.subvec(start, arma::size(center)) = center;
    }

   protected:
    virtual void encode(arma::dmat& genotype) = 0;
    void reset() noexcept { bed_.reset(); }

   private:
    BedReader bed_;
    double scale_factor_{};
    arma::rowvec center_;
};

class Grm : public IGrm
{
   public:
    explicit Grm(
        std::string_view bed_file,
        size_t chunk_size = DEFAULT_CHUNK_SIZE,
        const std::vector<std::string>& exclude_individuals = {});
    virtual arma::dmat compute();

   private:
    virtual arma::rowvec compute_center(const arma::dmat& genotype) = 0;
    void centerlize(arma::dmat& genotype);
    static double Scale(arma::dmat& grm);
};

class AddGrm : public Grm
{
    using Grm::Grm;

   private:
    void encode(arma::dmat& genotype) override;
    arma::rowvec compute_center(const arma::dmat& genotype) override;
};

class DomGrm : public Grm
{
    using Grm::Grm;

   private:
    void encode(arma::dmat& genotype) override;
    arma::rowvec compute_center(const arma::dmat& genotype) override;
};

}  // namespace gelex
