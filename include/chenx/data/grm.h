#pragma once

#include <armadillo>
#include <cstdint>
#include <string_view>
#include <vector>
#include "chenx/data/bed_reader.h"

namespace chenx
{

using arma::dmat;
using arma::rowvec;
void dom_encode(dmat& genotype);

class IGrm
{
   public:
    explicit IGrm(
        std::string_view bed_file,
        std::vector<std::string>&& exclude_individuals = {},
        uint64_t chunk_size = DEFAULT_CHUNK_SIZE);
    IGrm(const IGrm&) = delete;
    IGrm(IGrm&&) noexcept = default;
    IGrm& operator=(const IGrm&) = delete;
    IGrm& operator=(IGrm&&) noexcept = default;
    virtual ~IGrm() = default;

    double scale_factor() const noexcept { return scale_factor_; }
    void set_scale_factor(double scale_factor) { scale_factor_ = scale_factor; }

    const BedReader& bed() const noexcept { return bed_; }
    BedReader& bed() noexcept { return bed_; }

    const rowvec& center() const noexcept { return center_; }
    void set_center(rowvec&& center) { center_ = std::move(center); }
    void set_center(uint64_t start, const rowvec& center)
    {
        center_.subvec(start, arma::size(center)) = center;
    }

   protected:
    virtual void Encode(dmat& genotype) = 0;

   private:
    BedReader bed_;
    double scale_factor_{};
    rowvec center_;
};

class Grm : public IGrm
{
   public:
    explicit Grm(
        std::string_view bed_file,
        std::vector<std::string>&& exclude_individuals = {},
        uint64_t chunk_size = DEFAULT_CHUNK_SIZE);
    virtual dmat Compute();

   private:
    virtual rowvec ComputeCenter(const dmat& genotype) = 0;
    void Centerlize(dmat& genotype);
    static double Scale(dmat& grm);
};

class AddGrm : public Grm
{
    using Grm::Grm;

   private:
    void Encode(dmat& genotype) override;
    rowvec ComputeCenter(const dmat& genotype) override;
};

class DomGrm : public Grm
{
    using Grm::Grm;

   private:
    void Encode(dmat& genotype) override;
    rowvec ComputeCenter(const dmat& genotype) override;
};

}  // namespace chenx
